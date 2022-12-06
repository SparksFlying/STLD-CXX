#include "entity.h"
#include "utility.h"
#include <NTL/tools.h>

extern map<ERTreeEntry*, string> entryInfos;
extern TimerClock tc;
extern vector<double> costTimes;
extern vector<double> updateTimes;
extern vector<double> SLBCTimes;
extern vector<double> SBCTimes;
extern vector<double> SMETimes;
extern vector<double> SCGTimes;
extern double clearTime;

void test_buildRTree()
{
    PaillierFast crypto(1024);
    crypto.generate_keys();
    RTree t = deSeriRTree(config::dataFolder + "/" + "test_2_12_rtree.txt");
    shared_ptr<ERTreeEntry> root = encryptRTree(crypto, t);

    std::function<void(shared_ptr<ERTreeEntry>)> fn = [](shared_ptr<ERTreeEntry> entry){
        std::cout << entryInfos[entry.get()] << std::endl;
    };
    traverseERTree(root, fn);

    count(root, crypto);
}

void test_SIC()
{
    DSP server;
    assert(server.SIC(server.crypto->encrypt(1746), server.crypto->encrypt(113)).data == 0);  // 0
    assert(server.SIC(server.crypto->encrypt(1), server.crypto->encrypt(2)).data == 1); // 1
    assert(server.SIC(server.crypto->encrypt(2), server.crypto->encrypt(1)).data == 0);  // 0
    assert(server.SIC(server.crypto->encrypt(0), server.crypto->encrypt(0)).data == 1);  // 1
    assert(server.SIC(server.crypto->encrypt(1), server.crypto->encrypt(1000000)).data == 1);  // 1
    assert(server.SIC(server.crypto->encrypt(1000000), server.crypto->encrypt(1)).data == 0);  // 0
}

void test_SVC()
{
    DSP server;
    {
        Vec<Ciphertext> vec_a(NTL::INIT_SIZE_TYPE{}, 3, server.crypto->encrypt(1));
        Vec<Ciphertext> vec_b(NTL::INIT_SIZE_TYPE{}, 3, server.crypto->encrypt(1));
        string res = server.SVC(vec_a, vec_b, USE_PLAINTEXT)[0].data.data.get_str();
        assert(res[res.size() - 1] - 1 == '1');
        assert(res[res.size() - 2] - 1 == '1');
        assert(res[res.size() - 3] - 1 == '1');
    }
    {
        Vec<Ciphertext> vec_a(NTL::INIT_SIZE_TYPE{}, 3);
        Vec<Ciphertext> vec_b(NTL::INIT_SIZE_TYPE{}, 3);
        vec_a[0] = server.crypto->encrypt(2);
        vec_a[1] = server.crypto->encrypt(3);
        vec_a[2] = server.crypto->encrypt(4); 
        vec_b[0] = server.crypto->encrypt(3); 
        vec_b[1] = server.crypto->encrypt(1); 
        vec_b[2] = server.crypto->encrypt(2);
        string res = server.SVC(vec_a, vec_b, USE_PLAINTEXT)[0].data.data.get_str();
        assert(res[res.size() - 1] - 1 == '1');
        assert(res[res.size() - 2] - 1 == '0');
        assert(res[res.size() - 3] - 1 == '0');
    }
    {
        Vec<Ciphertext> vec_a(NTL::INIT_SIZE_TYPE{}, 32, server.crypto->encrypt(2));
        Vec<Ciphertext> vec_b(NTL::INIT_SIZE_TYPE{}, 32, server.crypto->encrypt(1));
        auto res = server.SVC(vec_a, vec_b, USE_PLAINTEXT);
        for(size_t i = 0; i < res.size() - 1; ++i){
            string str = res[i].data.data.get_str();
            for(char& c: str){
                assert(c - 1 == '0');
            }
        }
        size_t rest = size_t(vec_a.length()) % Vector::pack_count(32, *(server.crypto));
        string str = res.back().data.data.get_str();
        for(size_t j = 0; j < rest; ++j){
            assert(str[str.size() - 1 - j] - 1 == '0');
        }
    }
}

void test_SDDC()
{
    DSP server;
    {

    }
}

void test_getDistance()
{
    DSP server;
    EncRectType encRect(2, EncPointType(2));
    encRect[LEFT_BOTTOM_CORNER][0] = server.crypto->encrypt(1);
    encRect[LEFT_BOTTOM_CORNER][1] = server.crypto->encrypt(11);
    encRect[RIGHT_UP_CORNER][0] = server.crypto->encrypt(56);
    encRect[RIGHT_UP_CORNER][1] = server.crypto->encrypt(43);
    // min distance
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(0), server.crypto->encrypt(44)})) == 2);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(1), server.crypto->encrypt(43)})) == 0);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(28), server.crypto->encrypt(44)})) == 1);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(29), server.crypto->encrypt(43)})) == 0);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(57), server.crypto->encrypt(44)})) == 2);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(56), server.crypto->encrypt(43)})) == 0);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(57), server.crypto->encrypt(42)})) == 1);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(57), server.crypto->encrypt(10)})) == 2);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(28), server.crypto->encrypt(10)})) == 1);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(0), server.crypto->encrypt(0)})) == 1 + 11*11);
    assert(server.crypto->decrypt(server.getMinDistance(encRect, EncPointType{server.crypto->encrypt(0), server.crypto->encrypt(12)})) == 1);
    // max distance
    assert(server.crypto->decrypt(server.getMaxDistance(encRect, EncPointType{server.crypto->encrypt(28), server.crypto->encrypt(27)}, EncPointType{server.crypto->encrypt(0), server.crypto->encrypt(44)})) == 56*56 + 33*33);
    assert(server.crypto->decrypt(server.getMaxDistance(encRect, EncPointType{server.crypto->encrypt(28), server.crypto->encrypt(27)}, EncPointType{server.crypto->encrypt(0), server.crypto->encrypt(0)})) == 56*56 + 43*43);
    assert(server.crypto->decrypt(server.getMaxDistance(encRect, EncPointType{server.crypto->encrypt(28), server.crypto->encrypt(27)}, EncPointType{server.crypto->encrypt(57), server.crypto->encrypt(44)})) == 56*56 + 33*33);
    assert(server.crypto->decrypt(server.getMaxDistance(encRect, EncPointType{server.crypto->encrypt(28), server.crypto->encrypt(27)}, EncPointType{server.crypto->encrypt(57), server.crypto->encrypt(10)})) == 56*56 + 33*33);
}

void printResult(const vector<pair<vector<uint32_t>, int>>& W)
{
    for(auto& item: W){
        // string fmt = "point:[%s], score:%s";
        string fmt = format("point:[%s], score:%lu", point2str(item.first).c_str(), item.second);
        printf("%s\n", fmt.c_str());
    }
}

void test_STLD()
{
    DSP server;
    TimerClock tc;
    tc.update();
    server.loadRTree(config::dataFolder + "/" + "test_2_12_rtree.txt");
    EncPointType q(2);
    q[0] = server.crypto->encrypt(17);
    q[1] = server.crypto->encrypt(28);
    int k = 4;
    tc.update();
    auto W = server.topKQuery(q, k);
    printResult(W);
}

void experiment()
{
    DSP server;
    TimerClock tc;
    tc.update();
    server.loadRTree(config::dataFolder + "/" + "corr_2_1000_rtree.txt");
    printf("build time:%.3f\n", tc.getTimerMilliSec());

    size_t n_query = 10;
    double queryTime = 0;
    for(size_t i = 0; i < n_query; ++i){
        EncPointType q(2);
        q[0] = server.crypto->encrypt(Random::instance().rand_int(10000));
        q[1] = server.crypto->encrypt(Random::instance().rand_int(10000));
        int k = 2;
        tc.update();
        auto W = server.topKQuery(q, k);
        queryTime += tc.getTimerMilliSec();
    }
    printf("total query time:%.3f\n", queryTime);
    printf("cost times:");
    for(auto val: costTimes){
        std::cout << val << " ";
    }
    printf("\nSLBC times:");
    for(auto val: SLBCTimes){
        std::cout << val << " ";
    }
    printf("\nclear time:%.3f\n", clearTime);
}

void test_bitset()
{
    bitset<64> bs(1);
    assert(bs[0] == 1);
    assert(bs[63] == 0);
}

int main()
{
    // test_bitset();
    // experiment();
    // test_SIC();
    // test_SVC();
    // test_getDistance();
    test_STLD();
}