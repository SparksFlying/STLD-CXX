#include "entity.h"
#include "utility.h"

TimerClock tc;
vector<double> costTimes(4, 0);
vector<double> updateTimes(10, 0);
vector<double> SLBCTimes(2, 0);
vector<double> SBCTimes(7, 0);
vector<double> SMETimes(5, 0);
vector<double> SCGTimes(3, 0);
double clearTime = 0;


void DSP::recvKeys()
{
    rpc_client cli(config::CA_IP, config::CA_PORT);
    while(!cli.connect());
    size_t keySize = cli.call<size_t>("getKeySize");
    pair<string, string> param = cli.call<pair<string, string>>("getPub");
    array<string, 4> param2 = cli.call<array<string, 4>>("getPriv");
    PublicKey pk(keySize, Integer(param.first.c_str()), Integer(param.second.c_str()));
    PrivateKey sk(keySize, std::stoul(param2[0]), Integer(param2[1].c_str()), Integer(param2[2].c_str()), Integer(param2[3].c_str()));

    crypto = new PaillierFast(pk, sk);
    info("recv keys from %s:%d", config::CA_IP.c_str(), config::CA_PORT);
}

void DSP::loadRTree(const string& filePath)
{
    RTree t = deSeriRTree(filePath);
    root_ = encryptRTree(*crypto, t);
    count(root_, *crypto);
    
    // init dist cache
    extern map<ERTreeEntry*, string> entryInfos;
    for(auto& item: entryInfos){
        distCache[item.first] = vector<shared_ptr<Ciphertext>>(2, nullptr);
    }
}

void DSP::printEntryInfos()
{
    extern map<ERTreeEntry*, string> entryInfos;
    for(auto& item: entryInfos){
        string str = format(item.second.c_str(), crypto->decrypt(item.first->score_).get_str().c_str());
        printf("%s\n", str.c_str());
    }
}

vector<pair<vector<ValueType>, int>> DSP::topKQuery(const EncPointType& q, int k)
{
    vector<pair<EncPointType, Ciphertext>> W = STLD(root_, q, k);

    // returning
    vector<pair<vector<ValueType>, int>> ret(W.size());
    for(size_t i = 0; i < W.size(); ++i){
        ret[i].first.resize(W[i].first.size());
        for(size_t j = 0; j < W[i].first.size(); ++j){
            ret[i].first[j] = crypto->decrypt(W[i].first[j]).to_ulong();
        }
        ret[i].second = crypto->decrypt(W[i].second).to_ulong();
    }

    // clear
    tc.update();
    for(auto& item: distCache){
        for(size_t i = 0; i < 2; ++i){
            item.second[i].reset();
        }
    }
    traverseERTree(root_, this->clearFn);
    clearTime += tc.getTimerMilliSec();

    return ret;
}



vector<pair<EncPointType, Ciphertext>> DSP::SME(vector<pair<EncPointType, Ciphertext>>& X)
{
    if(X.size() < 2) return X;
    size_t k = X.size();
    vector<pair<size_t, string>> Y;
    vector<size_t> S;
    size_t bitLength = getMaxBitLength(k);
    vector<size_t> phi(k);
    for(size_t i = 0; i < k; ++i) phi[i] = i;
    Integer r1 = Random::instance().rand_int(1000) % Integer(k);
    for(size_t i = 0; i < k; ++i){
        S.push_back(phi[i]);
        Y.push_back({phi[i], (X[i].second.data * crypto->encrypt(r1).data).get_str()});
    }


    auto response = this->DAP_SME(std::move(Y));

    vector<vector<Ciphertext>> S_w(response.first.size(), vector<Ciphertext>(response.first[0].size()));
    vector<size_t> S_l = std::move(response.second);

    for(size_t i = 0; i < response.first.size(); ++i){
        for(size_t j = 0; j < response.first[i].size(); ++j){
            S_w[i][j].data = Integer(response.first[i][j].c_str());
        }
    }

    set<size_t> X_l;
    vector<pair<EncPointType, Ciphertext>> X_res(X.size());
    for(size_t j = 0; j < S_l.size(); ++j){
        for(size_t l = 0; l < k; ++l){
            if(S[l] == S_l[j]){
                X_res[j] = std::move(X[l]);
                X_l.insert(l);
            }
        }
    }

    vector<Ciphertext> P(2);
    vector<size_t> xs;
    {
        set<size_t> tmps(S.begin(), S.end());
        std::set_difference(tmps.begin(), tmps.end(), X_l.begin(), X_l.end(), std::inserter(xs, xs.begin()));
    }
    bitset<64> expr(S[xs[0]]);
    auto& expr2 = S_w[0];
    P[0] = expr2[0];

    if(expr[0] == 0){
        P[0] = SMinus(this->one_, P[0]);
    }
    for(size_t bidx = 1; bidx < bitLength; ++bidx){
        if(expr[bidx] == 0){
            P[0] = SM(P[0], SMinus(this->one_, expr2[bidx]));
        }else{
            P[0] = SM(P[0], expr2[bidx]);
        }
    }
    P[1] = SMinus(this->one_, P[0]);
    size_t dim = X[xs[0]].first.size();
    size_t x1 = xs[0], x2 = xs[1];
    X_res[S_l.size()].first.resize(dim);
    X_res[S_l.size() + 1].first.resize(dim);
    for(size_t didx = 0; didx < dim; ++didx){
        X_res[S_l.size()].first[didx] = SM(X[x1].first[didx], P[0]).data * SM(X[x2].first[didx], P[1]).data;
        X_res[S_l.size() + 1].first[didx] = SMinus(Ciphertext(X[x1].first[didx].data * X[x2].first[didx].data), X_res[S_l.size()].first[didx]);
    }
    X_res[S_l.size()].second = SM(X[x1].second, P[0]).data * SM(X[x2].second, P[1]).data;
    X_res[S_l.size() + 1].second = SMinus(Ciphertext(X[x1].second.data * X[x2].second.data), X_res[S_l.size()].second);
    return X_res;
}

Ciphertext DSP::SSED(const EncPointType& a, const EncPointType& b)
{
    Ciphertext dist = this->zero_;
    for(size_t didx = 0; didx < 2; ++didx){
        Ciphertext diff = SMinus(a[didx], b[didx]);
        dist.data = dist.data * SM(diff, diff).data;
    }
    return dist;
}

bool DSP::SDDC(const Ciphertext& da, const Ciphertext& db, const EncPointType& a, const EncPointType& b)
{
    // used to test
    // {
    //     Integer tmp_d_a = crypto->decrypt(da);
    //     Integer tmp_d_b = crypto->decrypt(db);
    //     bool domi = false;
    //     if(tmp_d_a > tmp_d_b) return false;
    //     else if(tmp_d_a < tmp_d_b) domi = true;
    //     for(size_t i = 2; i < a.size(); ++i){
    //         Integer tmp_a = crypto->decrypt(a[i]);
    //         Integer tmp_b = crypto->decrypt(b[i]);
    //         if(tmp_a > tmp_b) return false;
    //         else if(tmp_a < tmp_b) domi = true;
    //     }
    //     return domi;
    // }
    Ciphertext sum_a = da;
    Ciphertext sum_b = db;
    for(size_t i = 2; i < a.size(); ++i){
        sum_a.data = sum_a.data * a[i].data;
        sum_b.data = sum_b.data * b[i].data;
    }
    if(SIC(sum_b, sum_a, USE_PLAINTEXT).data == 1){
        return false;
    }
    Vec<Ciphertext> vec_a(NTL::INIT_SIZE_TYPE{}, a.size() - 1);
    Vec<Ciphertext> vec_b(NTL::INIT_SIZE_TYPE{}, b.size() - 1);
    vec_a[0] = da;
    vec_b[0] = db;
    for(size_t i = 2; i < a.size(); ++i){
        vec_a[i - 1] = a[i];
        vec_b[i - 1] = b[i];
    }
    auto cmpRes = SVC(vec_a, vec_b, USE_PLAINTEXT);
    for(size_t i = 0; i < cmpRes.size() - 1; ++i){
        string str = cmpRes[i].data.data.get_str();
        for(char& c: str){
            if(c  - 1 == '0') return false;
        }
    }
    size_t rest = size_t(vec_a.length()) % Vector::pack_count(32, *(this->crypto));
    string str = cmpRes.back().data.data.get_str();
    for(size_t j = 0; j < rest; ++j){
        if(str[str.size() - 1 - j] - 1 == '0') return false;
    }
    return true;
}

Ciphertext& DSP::getDistRelQ(shared_ptr<ERTreeEntry> obj, const EncPointType& q, bool isRect, bool minOrMax)
{
    if(isRect == true){
        if(!distCache[obj.get()][minOrMax]){
            Ciphertext dist;
            if(minOrMax == MIN){
                dist = getMinDistance(obj->rect_, q);
            }
            else{
                dist = getMaxDistance(obj->rect_, obj->midPoint_, q);
            }
            distCache[obj.get()][minOrMax].reset(new Ciphertext(std::move(dist)));
        }
        return *distCache[obj.get()][minOrMax];
    }
    else{
        if(!distCache[obj.get()][0]){
            distCache[obj.get()][0].reset(new Ciphertext(SSED(obj->data_, q)));   
        }
        return *distCache[obj.get()][0];
    }
}


bool _batch_count_leaf_check_in(DSP* dsp, const shared_ptr<ERTreeEntry>& p, const shared_ptr<ERTreeEntry>& o, const EncPointType& q)
{
    auto d_p = dsp->getDistRelQ(p, q);
    auto d_o_max = dsp->getDistRelQ(o, q, true, MAX);
    auto d_o_min = dsp->getDistRelQ(o, q, true, MIN);
    if(dsp->SDDC(d_p, d_o_max, p->data_, o->rect_[RIGHT_UP_CORNER]) == true &&
        dsp->SDDC(d_p, d_o_min, p->data_, o->rect_[LEFT_BOTTOM_CORNER]) == false){
        return true;
    }
    return false;
}

bool _batch_count_non_leaf_check_in(DSP* dsp, const shared_ptr<ERTreeEntry>& o_l, const shared_ptr<ERTreeEntry>& o, const EncPointType& q)
{
    auto d_o_l_min = dsp->getDistRelQ(o_l, q, true, MIN);
    auto d_o_max = dsp->getDistRelQ(o, q, true, MAX);
    auto d_o_l_max = dsp->getDistRelQ(o_l, q, true, MAX);
    auto d_o_min = dsp->getDistRelQ(o, q, true, MIN);
    if(dsp->SDDC(d_o_l_min, d_o_max, o_l->rect_[LEFT_BOTTOM_CORNER], o->rect_[RIGHT_UP_CORNER]) == true &&
        dsp->SDDC(d_o_l_max, d_o_min, o_l->rect_[RIGHT_UP_CORNER], o->rect_[LEFT_BOTTOM_CORNER]) == false)
    {
        return true;
    }
    return false;
}

bool DSP::_batch_count_check_in(vector<shared_ptr<ERTreeEntry>>& C, shared_ptr<ERTreeEntry> o, const EncPointType& q, bool isLeaf)
{
    if(isLeaf == true)
    {
        // seri
        for(auto& p: C){
            auto d_p = getDistRelQ(p, q);
            auto d_o_max = getDistRelQ(o, q, true, MAX);
            auto d_o_min = getDistRelQ(o, q, true, MIN);
            if(SDDC(d_p, d_o_max, p->data_, o->rect_[RIGHT_UP_CORNER]) == true &&
                SDDC(d_p, d_o_min, p->data_, o->rect_[LEFT_BOTTOM_CORNER]) == false){
                    return true;
            }
        }
        return false;
        // paral
        // vector<std::future<bool>> flags(C.size());
        // for(size_t i = 0; i < C.size(); ++i){
        //     flags[i] = std::async(_batch_count_leaf_check_in, this, C[i], o, q);
        // }
        // for(size_t i = 0; i < C.size(); ++i){
        //     bool flag = flags[i].get();
        //     if(flag) return true;
        // }
        // return false;
    }
    else
    {
        // seri
        for(auto& o_l: C){
            auto d_o_l_min = getDistRelQ(o_l, q, true, MIN);
            auto d_o_max = getDistRelQ(o, q, true, MAX);
            auto d_o_l_max = getDistRelQ(o_l, q, true, MAX);
            auto d_o_min = getDistRelQ(o, q, true, MIN);
            if(SDDC(d_o_l_min, d_o_max, o_l->rect_[LEFT_BOTTOM_CORNER], o->rect_[RIGHT_UP_CORNER]) == true &&
                SDDC(d_o_l_max, d_o_min, o_l->rect_[RIGHT_UP_CORNER], o->rect_[LEFT_BOTTOM_CORNER]) == false)
            {
                return true;
            }
        }
        return false;
        // paral
        // vector<std::future<bool>> flags(C.size());
        // for(size_t i = 0; i < C.size(); ++i){
        //     flags[i] = std::async(_batch_count_non_leaf_check_in, this, C[i], o, q);
        // }
        // for(size_t i = 0; i < C.size(); ++i){
        //     bool flag = flags[i].get();
        //     if(flag) return true;
        // }
        // return false;
    }
    return true;
}

void DSP::SLBC(const vector<shared_ptr<ERTreeEntry>>& Z, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C)
{
    for(auto& o: Z){
        tc.update();
        if(o->level_ > 2 && _batch_count_check_in(C, o, q, false)){
            SLBCTimes[0] += tc.getTimerMilliSec();
            SLBC(o->entries_, q, C);
        }else{
            SLBCTimes[0] += tc.getTimerMilliSec();
            tc.update();
            for(auto& o_l: C){
                auto d_o_l_min = getDistRelQ(o_l, q, true, MIN);
                auto d_o_max = getDistRelQ(o, q, true, MAX);
                if(SDDC(d_o_l_min, d_o_max, o_l->rect_[LEFT_BOTTOM_CORNER], o->rect_[RIGHT_UP_CORNER]) == true)
                {
                    o_l->score_.data = o_l->score_.data * o->count_;
                }
            }
            SLBCTimes[1] += tc.getTimerMilliSec();
        }
    }
}

void _SLBC_calc(DSP* dsp, const shared_ptr<ERTreeEntry>& o, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C)
{
    for(auto& o_l: C){
        auto d_o_l_min = dsp->getMinDistance(o_l->rect_, q);//this->getDistRelQ(o_l, q, true, MIN);
        auto d_o_max = dsp->getMaxDistance(o->rect_, o->midPoint_, q);//this->getDistRelQ(o, q, true, MAX);
        if(dsp->SDDC(d_o_l_min, d_o_max, o_l->rect_[LEFT_BOTTOM_CORNER], o->rect_[RIGHT_UP_CORNER]) == true)
        {
            o_l->score_.data = o_l->score_.data * o->count_;
        }
    }
}

void DSP::SLBCIter(const vector<shared_ptr<ERTreeEntry>>& Z, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C)
{
    vector<shared_ptr<ERTreeEntry>> curLevel = Z;
    vector<shared_ptr<ERTreeEntry>> nextLevel;
    vector<shared_ptr<ERTreeEntry>> needToCalc;
    tc.update();
    while(!curLevel.empty())
    {
        for(auto& o: curLevel){
            if(o->level_ > 2 && _batch_count_check_in(C, o, q, false)){
                for(auto& child: o->entries_){
                    nextLevel.push_back(child);
                }
            }
            else{
                needToCalc.push_back(o);
            }
        }
        curLevel = nextLevel;
        nextLevel.clear();
    }
    SLBCTimes[0] += tc.getTimerMilliSec();
    // printf("need to calc size:%lu\n", needToCalc.size());
    // seri
    tc.update();
    for(auto& o: needToCalc){
        for(auto& o_l: C){
            auto d_o_l_min = getDistRelQ(o_l, q, true, MIN);
            auto d_o_max = getDistRelQ(o, q, true, MAX);
            if(SDDC(d_o_l_min, d_o_max, o_l->rect_[LEFT_BOTTOM_CORNER], o->rect_[RIGHT_UP_CORNER]) == true)
            {
                // domiCache[o_l.get()][o.get()] = new bool(true);
                o_l->score_.data = o_l->score_.data * o->count_;
            }
            // else{
            //     domiCache[o_l.get()][o.get()] = new bool(false);
            // }
        }
    }
    SLBCTimes[1] += tc.getTimerMilliSec();
    // paral
    
    // vector<std::future<void>> tasks(needToCalc.size());
    // for(size_t i = 0; i < tasks.size(); ++i){
    //     tasks[i] = std::async(_SLBC_calc, this, needToCalc[i], q, std::ref(C));
    // }
    // for(size_t i = 0; i < tasks.size(); ++i){
    //     tasks[i].get();
    // }
}

void DSP::SBC(const vector<shared_ptr<ERTreeEntry>>& Z, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C)
{
    for(auto& o: Z){
        if(o->level_ > 1 && _batch_count_check_in(C, o, q, true)){
            if(o->level_ == 2){
                SBC(B(o), q, C);
            }else{
                SBC(o->entries_, q, C);
            }
        }else{
            if(o->level_ > 1){
                for(auto& p: C){
                    auto d_p = getDistRelQ(p, q);
                    auto d_o = getDistRelQ(o, q, true, MIN);
                    if(SDDC(d_p, d_o, p->data_, o->rect_[LEFT_BOTTOM_CORNER]) == true){
                        p->score_.data = p->score_.data * o->count_;
                    }
                }
            }else{
                for(auto& p: C){
                    auto d_p = getDistRelQ(p, q);
                    auto d_o = getDistRelQ(o, q);
                    if(SDDC(d_p, d_o, p->data_, o->data_) == true){
                        p->score_.data = p->score_.data * o->count_;
                    }
                }
            }
        }
    }
}

bool DSP::_topk_check_in(const vector<EncPointType>& F, const EncPointType& q, shared_ptr<ERTreeEntry> o, bool type)
{
    for(auto& p: F){
        auto d_p = SSED(p, q);
        auto d_o = getDistRelQ(o, q, true, MIN);
        if(SDDC(d_p, d_o, p, o->rect_[LEFT_BOTTOM_CORNER]) == true){
            return type == EXIST;
        }
    }
    return type != EXIST;
}

bool cmpFunctor(shared_ptr<ERTreeEntry>& a, shared_ptr<ERTreeEntry>& b)
{
    return *a < *b;
}

vector<pair<EncPointType, Ciphertext>> DSP::STLD(shared_ptr<ERTreeEntry> root, const EncPointType& q, int k)
{
    priority_queue<shared_ptr<ERTreeEntry>, vector<shared_ptr<ERTreeEntry>>, decltype(&cmpFunctor)> H(cmpFunctor);
    vector<EncPointType> F, F_l;
    vector<shared_ptr<ERTreeEntry>> T;
    Ciphertext phi = this->zero_;
    vector<shared_ptr<ERTreeEntry>> roots{root_};
    tc.update();
    SLBC(roots, q, roots);
    costTimes[2] += tc.getTimerMilliSec();
    vector<pair<EncPointType, Ciphertext>> W(k);
    size_t dim = root_->rect_[LEFT_BOTTOM_CORNER].size();
    for(size_t i = 0; i < k; ++i){
        W[i].first.resize(dim);
        for(size_t j = 0; j < dim; ++j){
            W[i].first[j].data = this->const_max_.data;
        }
        W[i].second = crypto->encrypt(Integer(-1));
    }

    H.push(root);

    while(!H.empty()){
        auto o = H.top();
        H.pop();
        tc.update();
        if(SIC(o->score_, phi, USE_PLAINTEXT).data == 1 && _topk_check_in(F, q, o, EXIST)){
            break;
        }
        costTimes[0] += tc.getTimerMilliSec();
        if(o->level_ > 2){
            auto& Z = o->entries_;
            vector<shared_ptr<ERTreeEntry>> C;
            tc.update();
            for(auto& o_c: Z){
                if(_topk_check_in(F, q, o_c, ANY)){
                    C.push_back(o_c);
                }
                H.push(o_c);
            }
            costTimes[1] += tc.getTimerMilliSec();
            tc.update();
            SLBC(roots, q, C);
            costTimes[2] += tc.getTimerMilliSec();
        }else{
            auto& Z = B(o);
            tc.update();
            UpdateResult(W, F, F_l, Z, q, phi, root_, k, T, threshold_);
            costTimes[3] += tc.getTimerMilliSec();
            // printEntryInfos();
        }
    }
    tc.update();
    UpdateResult(W, F, F_l, {}, q, phi, root_, k, T, 0);
    costTimes[3] += tc.getTimerMilliSec();
    // printEntryInfos();
    return W;
}

bool DSP::_update_check_in(const vector<EncPointType>& F, const EncPointType& q, const EncPointType& p_e_data, bool type)
{
    for(auto& p: F){
        auto d_p = SSED(p, q);
        auto d_p_e = SSED(p_e_data, q);
        if(SDDC(d_p, d_p_e, p, p_e_data) == true){
            return  type == EXIST;
        }
    }
    return type != EXIST;
}

void DSP::UpdateResult(vector<pair<EncPointType, Ciphertext>>& W,
                        vector<EncPointType>& F,
                        vector<EncPointType>& F_l,
                        const vector<shared_ptr<ERTreeEntry>>& Z,
                        const EncPointType& q,
                        Ciphertext& phi,
                        shared_ptr<ERTreeEntry> root,
                        int k,
                        vector<shared_ptr<ERTreeEntry>>& T,
                        int delta_0)
{
    // seri
    for(auto& p_e: Z){
        if(_update_check_in(F, q, p_e->data_, ANY) == 1){
            T.push_back(p_e);
        }
    }
    // paral
    // {
    //     vector<std::future<bool>> flags(Z.size());
    //     for(size_t i = 0; i < Z.size(); ++i){
    //         flags[i] = std::async(&DSP::_update_check_in, this, F, q, Z[i]->data_, ANY);
    //     }
    //     for(size_t i = 0; i < Z.size(); ++i){
    //         bool flag = flags[i].get();
    //         if(flag) T.push_back(Z[i]);
    //     }
    // }
    if(T.size() > delta_0){
        SBC(vector<shared_ptr<ERTreeEntry>>{root}, q, T);
        for(auto& p_e: T){
            Ciphertext epsilon = SIC(p_e->score_, phi, USE_CIPHERTEXT);
            EncPointType o_l(p_e->data_.size());
            for(size_t i = 0; i < o_l.size(); ++i){
                o_l[i].data = SM(epsilon, p_e->data_[i]).data * SM(SMinus(this->one_, epsilon), W[k - 1].first[i]).data;
                W[k - 1].first[i].data = SMinus(Ciphertext(p_e->data_[i].data * W[k - 1].first[i].data), o_l[i]).data;
            }
            Ciphertext o_l_score = SM(epsilon, p_e->score_).data * SM(SMinus(this->one_, epsilon), W[k - 1].second).data;
            W[k - 1].second.data = SMinus(Ciphertext(p_e->score_.data * W[k - 1].second.data), o_l_score).data;

            if(SIC(o_l_score, this->zero_, USE_PLAINTEXT).data == 0 && _update_check_in(F, q, o_l, ANY)){
                F_l.push_back(o_l);
            }
            W = SME(W);
            phi.data = W[k - 1].second.data;
            F = F_l;
            F.push_back(W[k - 1].first);
        }
        T.clear();
    }
}


Ciphertext DSP::getMinDistance(const EncRectType& encRect, const EncPointType& q)
{
    Ciphertext dist = this->zero_;
    vector<bool> lowFlags(2, false);
    vector<bool> highFlags(2, false);

    for(size_t i = 0; i < 2; ++i){
        lowFlags[i] = SIC(encRect[LEFT_BOTTOM_CORNER][i], q[i], USE_PLAINTEXT).data.to_int();
        highFlags[i] = SIC(q[i], encRect[RIGHT_UP_CORNER][i], USE_PLAINTEXT).data.to_int();
    }

    if(!lowFlags[0] && !lowFlags[1]){
        dist.data = SSED(encRect[LEFT_BOTTOM_CORNER], q).data;
    }
    else if(!lowFlags[0] && lowFlags[1] && highFlags[1]){
        Ciphertext df0 = SMinus(encRect[LEFT_BOTTOM_CORNER][0], q[0]);
        dist.data = SM(df0, df0).data;
    }
    else if(!lowFlags[0] && !highFlags[1]){
        Ciphertext df0 = SMinus(encRect[LEFT_BOTTOM_CORNER][0], q[0]);
        Ciphertext df1 = SMinus(encRect[RIGHT_UP_CORNER][1], q[1]);
        dist.data = SM(df0, df0).data * SM(df1, df1).data;
    }
    else if(lowFlags[0] && highFlags[0] && !highFlags[1]){
        Ciphertext df1 = SMinus(q[1], encRect[RIGHT_UP_CORNER][1]);
        dist.data = SM(df1, df1).data;
    }
    else if(!highFlags[0] && !highFlags[1]){
        dist.data = SSED(encRect[RIGHT_UP_CORNER], q).data;
    }
    else if(!highFlags[0] && lowFlags[1] && highFlags[1]){
        Ciphertext df0 = SMinus(q[0], encRect[RIGHT_UP_CORNER][0]);
        dist.data = SM(df0, df0).data;
    }
    else if(!highFlags[0] && !lowFlags[1]){
        Ciphertext df0 = SMinus(encRect[RIGHT_UP_CORNER][0], q[0]);
        Ciphertext df1 = SMinus(encRect[LEFT_BOTTOM_CORNER][1], q[1]);
        dist.data = SM(df0, df0).data * SM(df1, df1).data;
    }
    else if(lowFlags[0] && highFlags[0] && !lowFlags[1]){
        Ciphertext df1 = SMinus(encRect[LEFT_BOTTOM_CORNER][1], q[1]);
        dist.data = SM(df1, df1).data;
    }
    return dist;
}

Ciphertext DSP::getMaxDistance(const EncRectType& encRect, const EncPointType& midPoint, const EncPointType& q)
{
    bool xFlag = SIC(q[0], midPoint[0], USE_PLAINTEXT).data.to_int();
    bool yFlag = SIC(q[1], midPoint[1], USE_PLAINTEXT).data.to_int();
    if(!xFlag && !yFlag){
        return SSED(q, encRect[LEFT_BOTTOM_CORNER]);
    }
    else if(!xFlag && yFlag){
        return SSED(q, EncPointType{encRect[LEFT_BOTTOM_CORNER][0], encRect[RIGHT_UP_CORNER][1]});
    }
    else if(xFlag && !yFlag){
        return SSED(q, EncPointType{encRect[RIGHT_UP_CORNER][0], encRect[LEFT_BOTTOM_CORNER][1]});
    }
    return SSED(q, encRect[RIGHT_UP_CORNER]);
}


// Secure Integer Comparison Protocol
// if a <= b then return 1, else return 0
Ciphertext DSP::SIC(const Ciphertext& a, const Ciphertext& b, bool level) {
    Ciphertext X = _times(a, 2);
    Ciphertext Y = _times(b, 2).data * this->one_.data;
    Integer coin = Random::instance().rand_int(crypto->get_pub().n) % Integer(2);
    Ciphertext Z;
    if(coin == 1){
        Z = SMinus(X, Y);
    }else{
        Z = SMinus(Y, X);
    }

    int maxBitLength = getMaxBitLength(crypto->get_pub().n) / 4 - 2;
    Integer r = Random::instance().rand_int(crypto->get_pub().n - 1) % (Integer(2).pow(maxBitLength)) + 1;
    Ciphertext c = _times(Z, r);

    // Ciphertext ret = Integer(dap_.call<string>("SIC", c.data.get_str()).c_str());
    string response = this->DAP_SIC(c.data.get_str(), level);
    if(level == USE_PLAINTEXT){
        Ciphertext ret;
        ret.data = Integer(response.c_str());
        if(coin == 0){
            ret.data = 1 - ret.data;
        }
        return ret;
    }
    Ciphertext ret = Integer(response.c_str());
    if(coin == 0){
        ret = SMinus(this->one_, ret);
    }
    return ret;
}

// packed SIC
vector<PackedCiphertext> DSP::SVC(const Vec<Ciphertext>& a, const Vec<Ciphertext>& b, bool level){
    //debug("dsp SVC start");
    assert(a.length() == b.length());
    size_t n = a.length();
    size_t numIntegerPerCiphertext = Vector::pack_count(32, *(this->crypto));
    size_t numPackedCiphertext = size_t(std::ceil(double(n) / double(numIntegerPerCiphertext)));
    
    Integer coin = Random::instance().rand_int(crypto->get_pub().n) % Integer(2);
    Vec<Ciphertext> Z(NTL::INIT_SIZE_TYPE{}, n);

    if(coin == 1){
        for(size_t i = 0; i < n; ++i){
            Z[i] = SMinus(_times(a[i], 2), _times(b[i], 2).data * this->one_.data).data * crypto->encrypt(Integer(2).pow(32 - 1)).data;
        }
    }else{
        for(size_t i = 0; i < n; ++i){
            Z[i] = SMinus(_times(b[i], 2).data * this->one_.data, _times(a[i], 2)).data * crypto->encrypt(Integer(2).pow(32 - 1)).data;
        }
    }

    Integer r = Random::instance().rand_int_bits((getMaxBitLength(crypto->get_pub().n) / 4) - 2);

    r = 1;
    for(auto& z : Z){
        z = _times(z, r);
    }

    vector<PackedCiphertext> C(numPackedCiphertext, PackedCiphertext(this->zero_, numIntegerPerCiphertext, 32));

    size_t idx = 0;
    for(size_t i = 0; i < numPackedCiphertext; ++i){
        size_t s = idx;
        size_t e = std::min(idx + numIntegerPerCiphertext, n);
        Vec<Ciphertext> packedC(NTL::INIT_SIZE_TYPE{}, numIntegerPerCiphertext, this->zero_);
        for(size_t j = s; j < e; ++j, ++idx){
            packedC[j - s].data = Z[idx].data;
        }
        C[i] = Vector::pack_ciphertexts(packedC, 32, *crypto);
    }

    vector<string> seriC(C.size());
    for(size_t i = 0; i < C.size(); ++i){
        seriC[i] = C[i].data.data.get_str();
    }
    // vector<string> packedCmpResult = dap_.call<vector<string>>("SVC", seriC);
    vector<string> packedCmpResult = this->DAP_SVC(seriC, level);
    if(level == USE_PLAINTEXT){
        vector<PackedCiphertext> ret(numPackedCiphertext);
        for(size_t i = 0; i < numPackedCiphertext; ++i){
            if(coin == 0){
                for(size_t j = 0; j < packedCmpResult[i].size(); ++j){
                    packedCmpResult[i][j] = '0' + ('1' - (packedCmpResult[i][j] - 1)) + 1;
                }
            }
            ret[i].data.data = Integer(packedCmpResult[i].c_str());
        }
        return ret;
    }

    vector<PackedCiphertext> ret(numPackedCiphertext);
    for(size_t i = 0; i < numPackedCiphertext; ++i){
        if(coin == 1){
            ret[i] = PackedCiphertext(this->zero_, numIntegerPerCiphertext, 32);
            ret[i].data.data = Integer(packedCmpResult[i].c_str());
        }else{
            ret[i] = Vector::encrypt_pack(Vec<Integer>(NTL::INIT_SIZE_TYPE{}, numIntegerPerCiphertext, 1), 32, *crypto);
            ret[i].data.data = SMinus(ret[i].data.data, Integer(packedCmpResult[i].c_str())).data;
        }
    }
    return ret;
}

string DSP::DAP_SIC(const string& c, bool level){
    //debug("call SIC from %s", conn.lock()->remote_address().c_str());
    Integer m = crypto->decrypt(Ciphertext(Integer(c.c_str()))) % crypto->get_pub().n;
    Integer u = 0;
    if(getMaxBitLength(m >= 0 ? m : -m) > getMaxBitLength(crypto->get_pub().n) / 2){
        u = 1;
    }
    if(level == USE_PLAINTEXT){
        return u.get_str();
    }else{
        return crypto->encrypt(u).data.get_str();
    }
}

vector<string> DSP::DAP_SVC(const vector<string>& C, bool level){
    //debug("call SVC from %s", conn.lock()->remote_address().c_str());
    vector<Vec<Integer>> M(C.size());
    vector<string> ret(M.size());
    size_t numIntegerPerCiphertext = Vector::pack_count(32, *(this->crypto));
    for(size_t i = 0; i < C.size(); ++i){
        PackedCiphertext packed(Integer(C[i].c_str()), numIntegerPerCiphertext, 32);
        M[i] = Vector::decrypt_pack(packed, *(this->crypto));

        Vec<Integer> packedCmpResult(NTL::INIT_SIZE_TYPE{}, M[i].length(), 0);
        for(size_t j = 0; j < numIntegerPerCiphertext; ++j){
            M[i][j] -= Integer(2).pow(32 - 1);
            if(getMaxBitLength(M[i][j] % crypto->get_pub().n) > getMaxBitLength(crypto->get_pub().n) / 2){
                packedCmpResult[j] = 1;
            }
        }
        if(level == USE_CIPHERTEXT){
            ret[i] = PackedCiphertext(Vector::encrypt_pack(packedCmpResult, 32, *crypto)).data.data.get_str();
        }else{
            for(size_t k = 0; k < M[i].length(); ++k){
                ret[i]  = string(1, packedCmpResult[k].get_str()[0] + 1) + ret[i];
            }
        }
    }
    return ret;
}

// secure multiply protocol
// res = E(a*b)
Ciphertext DSP::SM(const Ciphertext& a, const Ciphertext& b){
    Random& r = Random::instance();
    Integer ra = r.rand_int(crypto->get_pub().n);
    Integer rb = r.rand_int(crypto->get_pub().n);
    Ciphertext tmp_a = crypto->encrypt(ra).data * a.data;
    Ciphertext tmp_b = crypto->encrypt(rb).data * b.data;
    // Ciphertext h = Integer(dap_.call<string>("SM", tmp_a.data.get_str(), tmp_b.data.get_str()).c_str());
    Ciphertext h = Integer(this->DAP_SM(tmp_a.data.get_str(), tmp_b.data.get_str()).c_str());

    auto s = h.data * a.data.pow_mod_n(crypto->get_pub().n - rb, *crypto->get_n2());
    auto tmp_s = s * b.data.pow_mod_n(crypto->get_pub().n - ra, *crypto->get_n2());
    Ciphertext res = tmp_s * crypto->encrypt(ra * rb % crypto->get_pub().n).data.pow_mod_n(crypto->get_pub().n - 1, *crypto->get_n2());
    return res;
}

string DSP::DAP_SM(const string& a, const string& b){
    //debug("call SM from %s", conn.lock()->remote_address().c_str());
    //printf("call SM from %s", conn.lock()->remote_address().c_str());
    Integer ha = crypto->decrypt(Integer(a.c_str()) % *crypto->get_n2());
    Integer hb = crypto->decrypt(Integer(b.c_str()) % *crypto->get_n2());
    Integer h = (ha * hb) % crypto->get_pub().n;
    return crypto->encrypt(h).data.get_str();
}

pair<vector<vector<string>>, vector<size_t>> DSP::DAP_SME(vector<pair<size_t, string>>&& Y)
{
    size_t n = Y.size();
    vector<pair<size_t, Integer>> X(n);
    for(size_t i = 0; i < n; ++i){
        X[i].first = Y[i].first;
        X[i].second = crypto->decrypt(Ciphertext(Integer(Y[i].second.c_str())));
    }
    Integer bitLength = getMaxBitLength(n);
    sort(X.begin(), X.end(), [](const pair<size_t, Integer>& a, const pair<size_t, Integer>& b){
        return a.second > b.second;
    });
    Integer r2 = Random::instance().rand_int(n) % Integer(n - 1);
    array<Integer, 2> selectedIdx{r2, n - 1};
    vector<vector<string>> S_w(2, vector<string>(bitLength.to_ulong()));
    for(size_t i = 0; i < 2; ++i){
        bitset<64> bs(X[selectedIdx[i].to_ulong()].first);
        for(size_t j = 0; j < bitLength.to_ulong(); ++j){
            int bit = bs[j];
            S_w[i][j] = crypto->encrypt(Integer(bit)).data.get_str();
        }
    }
    vector<size_t> S_l(n - 2);
    for(size_t i = 0, idx = 0; i < n; ++i){
        if(i != selectedIdx[0] && i != selectedIdx[1]){
            S_l[idx++] = X[i].first;
        }
    }
    std::random_shuffle(S_l.begin(), S_l.end());
    return pair<vector<vector<string>>, vector<size_t>>{S_w, S_l};
}