#pragma once
#include "config.hpp"

#define LEFT_BOTTOM_CORNER 0
#define RIGHT_UP_CORNER 1


using RectType = vector<vector<Integer>>;
using EncRectType = vector<vector<Ciphertext>>;
using PointType = vector<Integer>;
using EncPointType = vector<Ciphertext>;
using ValueType = uint32_t;

class RTreeEntry
{
public:
    RTreeEntry(){}
public:
    bool is_leaf;
    RectType rect;
    PointType data;
    vector<shared_ptr<RTreeEntry>> entries;
};
class RTree
{
public:
    RTree(){}
    inline int getHeight(){
        return this->height;
    }
public:
    shared_ptr<RTreeEntry> root;
    int height;
};

class ERTreeEntry
{
public:
    static std::function<bool(const Ciphertext& a, const Ciphertext& b)> cmpFunctor;
    ERTreeEntry(const PaillierFast& cryoto, const RTreeEntry& entry);
    inline bool operator < (const ERTreeEntry& other){
        return !ERTreeEntry::cmpFunctor(this->score_, other.score_);
    }
public:
    int level_;
    bool is_leaf_;
    vector<shared_ptr<ERTreeEntry>> entries_;
    EncRectType rect_;
    EncPointType midPoint_;
    EncPointType data_;
    Integer count_; // actually as ciphertext
    Ciphertext score_;
};

extern size_t getMaxBitLength(const Integer& a);

extern void splitString(const std::string& s, std::vector<std::string>& v, const std::string& c);

extern string lstrip(const string & str, const string & chars = " ");

extern Integer min(const Integer& a, const Integer& b);

extern Integer max(const Integer& a, const Integer& b);

extern EncRectType encryptRectanle(const PaillierFast& crypto, const RectType& rect);

extern RectType decryptRectangle(const PaillierFast& crypto, const EncRectType& encRect);

extern RectType getMBR(const vector<RectType>& rects);

extern RTree deSeriRTree(const string& filePath);

extern shared_ptr<ERTreeEntry> encryptRTree(const PaillierFast& crypto, RTree& t);

extern void traverseERTree(shared_ptr<ERTreeEntry> root, std::function<void(shared_ptr<ERTreeEntry>)>& fn);

extern void count(shared_ptr<ERTreeEntry> root, const PaillierFast& crypto);

extern vector<PointType> read(const string& path);

extern vector<shared_ptr<ERTreeEntry>>& B(shared_ptr<ERTreeEntry> entry);

template<class... T>
string format(const char *fmt, const T&...t)
{
    const auto len = snprintf(nullptr, 0, fmt, t...);
    string r;
    r.resize(static_cast<size_t>(len) + 1);
    snprintf(&r.front(), len + 1, fmt, t...);  // Bad boy
    r.resize(static_cast<size_t>(len));
 
    return r;
}

template<class T>
string point2str(const vector<T>& p)
{
    string content;
    for(size_t i = 0; i < p.size() - 1; ++i){
        content += Integer(p[i]).get_str() + ",";
    }
    content += Integer(p.back()).get_str();
    return content;
}

/*===================================================================================================*/
/*============================================== Timer ==============================================*/
/*===================================================================================================*/
// provides second, millisecond and microsecond timing capability
class TimerClock
{
public:
    TimerClock()
    {
        update();
    }
    
    ~TimerClock(){}
    
    void update()
    {
        _start = high_resolution_clock::now();
    }
    double getTimerSecond()
    {
        return getTimerMicroSec() * 0.000001;
    }
    double getTimerMilliSec()
    {
        return getTimerMicroSec()*0.001;
    }
    long long getTimerMicroSec()
    {
        return duration_cast<microseconds>(high_resolution_clock::now() - _start).count();
    }
private:
    time_point<high_resolution_clock>_start;
};