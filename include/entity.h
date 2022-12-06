#pragma once
#include "utility.h"
#include "logging.h"

#define USE_PLAINTEXT 0
#define USE_CIPHERTEXT 1
#define EXIST 1
#define ANY 0
#undef MIN
#undef MAX
#define MIN 0
#define MAX 1

inline void _scores_clear(shared_ptr<ERTreeEntry> entry, const Ciphertext& val)
{
    entry->score_.data = val.data;
}

class CA : public rpc_server{
public:
    CA(const uint32_t key_size = config::KEY_SIZE) : rpc_server(config::CA_PORT, 1), crypto(new PaillierFast(key_size)){
        crypto->generate_keys();
        this->register_handler("getKeySize", &CA::getKeySize, this);
        this->register_handler("getPub", &CA::getPub, this);
        this->register_handler("getPriv", &CA::getPriv, this);
    }
    void saveKeys(const string& filePath){
        std::ofstream out(filePath, std::ios::out | std::ios::binary);
        if(!out.is_open()){
            // info("open %s failed!", filePath.c_str());
            return;
        }
        // save key size
        out << crypto->get_pub().key_size_bits << ' ';
        // save pub key
        out << crypto->get_pub().n.get_str() << ' ' << crypto->get_pub().g.get_str() << ' ';
        // save priv key
        out << std::to_string(crypto->get_priv().a_bits) << ' ';
        out << crypto->get_priv().p.get_str() << ' ';
        out << crypto->get_priv().q.get_str() << ' ';
        out << crypto->get_priv().a.get_str();
        out.close();
        // info("save done\n");
    }
    CA(const string& filePath) : rpc_server(config::CA_PORT, 1), crypto(nullptr){
        readKeys(filePath);
        this->register_handler("getKeySize", &CA::getKeySize, this);
        this->register_handler("getPub", &CA::getPub, this);
        this->register_handler("getPriv", &CA::getPriv, this);
    }
    void readKeys(const string& filePath){
        std::ifstream in(filePath, std::ios::in | std::ios::binary);
        if(!in.is_open()){
            // info("open %s failed!", filePath.c_str());
            return;
        }
        size_t keySize;
        pair<string, string> param;
        array<string, 4> param2;

        in >> keySize;
        in >> param.first >> param.second;
        for(size_t i = 0; i < 4; ++i){
            in >> param2[i];
        }
        PublicKey pk(keySize, Integer(param.first.c_str()), Integer(param.second.c_str()));
        PrivateKey sk(keySize, std::stoul(param2[0]), Integer(param2[1].c_str()), Integer(param2[2].c_str()), Integer(param2[3].c_str()));

        if(crypto) delete crypto;
        crypto = new PaillierFast(pk, sk);
        // info("read done\n");
    }
    size_t getKeySize(rpc_conn conn){
        return crypto->get_pub().key_size_bits;
    }
    pair<string, string> getPub(rpc_conn conn){
        return {crypto->get_pub().n.get_str(), crypto->get_pub().g.get_str()};
    }
    array<string, 4> getPriv(rpc_conn conn){
        return {std::to_string(crypto->get_priv().a_bits), crypto->get_priv().p.get_str(), crypto->get_priv().q.get_str(), crypto->get_priv().a.get_str()};
    }

private:
    PaillierFast* crypto;
};

class DSP : public rpc_server{
public:
    DSP(int threshold = 4) : rpc_server(config::DSP_PORT, config::NUM_THREADS), crypto(nullptr), dap_(config::DAP_IP, config::DAP_PORT), threshold_(threshold){
        recvKeys();

        this->zero_ = crypto->encrypt(0);
        this->one_ = crypto->encrypt(1);
        this->const_max_ = crypto->encrypt(Integer(2).pow(64));

        clearFn = std::bind(_scores_clear, std::placeholders::_1, this->zero_);

        dap_.enable_auto_reconnect(true);
        dap_.enable_auto_heartbeat(true);
        // while(!dap_.connect());

        // this->register_handler("topKQuery", &DSP::topKQuery, this);
        decltype(ERTreeEntry::cmpFunctor) cmpFunctor = [this](const Ciphertext& a, const Ciphertext& b) -> bool{
            if(this->SIC(a, b, USE_PLAINTEXT).data.get_str()[0] == '1'){
                return true;
            }
            return false;
        };
        ERTreeEntry::cmpFunctor = cmpFunctor;
    }

    void recvKeys();
    vector<pair<vector<uint32_t>, int>> topKQuery(const EncPointType& q, int k);
    
    // DSP protocols
    vector<pair<EncPointType, Ciphertext>> SME(vector<pair<EncPointType, Ciphertext>>& X);
    Ciphertext SIC(const Ciphertext& a, const Ciphertext& b, bool level = USE_PLAINTEXT);
    vector<PackedCiphertext> SVC(const Vec<Ciphertext>& a, const Vec<Ciphertext>& b, bool level = USE_PLAINTEXT);
    Ciphertext SM(const Ciphertext& a, const Ciphertext& b);
    Ciphertext SSED(const EncPointType& a, const EncPointType& b);
    bool SDDC(const Ciphertext& da, const Ciphertext& db, const EncPointType& a, const EncPointType& b);
    Ciphertext& getDistRelQ(shared_ptr<ERTreeEntry> obj, const EncPointType& q, bool isRect = false, bool minOrMax = MIN);
    bool _batch_count_check_in(vector<shared_ptr<ERTreeEntry>>& C, shared_ptr<ERTreeEntry> o, const EncPointType& q, bool isLeaf = false);
    bool _topk_check_in(const vector<EncPointType>& F, const EncPointType& q, shared_ptr<ERTreeEntry> o, bool type = EXIST);
    bool _update_check_in(const vector<EncPointType>& F, const EncPointType& q, const EncPointType& p_e_data, bool type = ANY);
    void SLBC(const vector<shared_ptr<ERTreeEntry>>& Z, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C);
    void SLBCIter(const vector<shared_ptr<ERTreeEntry>>& Z, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C);
    void SBC(const vector<shared_ptr<ERTreeEntry>>& Z, const EncPointType& q, vector<shared_ptr<ERTreeEntry>>& C);
    vector<pair<EncPointType, Ciphertext>> STLD(shared_ptr<ERTreeEntry> root, const EncPointType& q, int k);
    void UpdateResult(vector<pair<EncPointType, Ciphertext>>& W,
                        vector<EncPointType>& F,
                        vector<EncPointType>& F_l,
                        const vector<shared_ptr<ERTreeEntry>>& Z,
                        const EncPointType& q,
                        Ciphertext& phi,
                        shared_ptr<ERTreeEntry> root,
                        int k,
                        vector<shared_ptr<ERTreeEntry>>& T,
                        int delta_0);
    Ciphertext getMinDistance(const EncRectType& encRect, const EncPointType& q);
    Ciphertext getMaxDistance(const EncRectType& encRect, const EncPointType& midPoint, const EncPointType& q);
    void loadRTree(const string& filePath);
    void printEntryInfos();

    inline Ciphertext _times(const Ciphertext& a, Integer times){
        return a.data.pow_mod_n(times, *(crypto->get_n2()));
    }

    // res=E(a-b)
    inline Ciphertext SMinus(const Ciphertext& a,const Ciphertext& b) {
        return a.data * b.data.pow_mod_n(crypto->get_pub().n - 1, *(crypto->get_n2()));
    }

    // DAP protocols
    pair<vector<vector<string>>, vector<size_t>> DAP_SME(vector<pair<size_t, string>>&& Y);
    string DAP_SIC(const string& c, bool level);
    vector<string> DAP_SVC(const vector<string>& C, bool level);
    string DAP_SM(const string& a, const string& b);

public:
    PaillierFast* crypto;
    rpc_client dap_;
    shared_ptr<ERTreeEntry> root_;
    int threshold_ = 4;
    map<ERTreeEntry*, vector<shared_ptr<Ciphertext>>> distCache;
    // <ERTreeEntry*, map<ERTreeEntry*, bool*>> domiCache;

    // constant
    Ciphertext zero_;
    Ciphertext one_;
    Ciphertext const_max_;

    // clear func
    std::function<void(shared_ptr<ERTreeEntry>)> clearFn;
};
