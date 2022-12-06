#include "utility.h"
#include <fstream>

#define LEFTSTRIP 0
#define RIGHTSTRIP 1
#define BOTHSTRIP 2

map<ERTreeEntry*, string> entryInfos;

std::function<bool(const Ciphertext& a, const Ciphertext& b)> ERTreeEntry::cmpFunctor;

void splitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while(std::string::npos != pos2)
	{
        auto sub = s.substr(pos1, pos2-pos1);
		if(!sub.empty()) v.push_back(sub);
	
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if(pos1 != s.length())
    {
        auto sub = s.substr(pos1);
        if(!sub.empty()) v.push_back(s.substr(pos1));
    }
}

std::string do_strip(const std::string &str, int striptype, const std::string&chars)
{
    std::string::size_type strlen = str.size();
    std::string::size_type charslen = chars.size();
    std::string::size_type i, j;

    //默认情况下，去除空白符
    if (0 == charslen)
    {
        i = 0;
        //去掉左边空白字符
        if (striptype != RIGHTSTRIP)
        {
            while (i < strlen&&::isspace(str[i]))
            {
                i++;
            }
        }
        j = strlen;
        //去掉右边空白字符
        if (striptype != LEFTSTRIP)
        {
            j--;
            while (j >= i&&::isspace(str[j]))
            {
                j--;
            }
            j++;
        }
    }
    else
    {
        //把删除序列转为c字符串
        const char*sep = chars.c_str();
        i = 0;
        if (striptype != RIGHTSTRIP)
        {
            //memchr函数：从sep指向的内存区域的前charslen个字节查找str[i]
            while (i < strlen&&memchr(sep, str[i], charslen))
            {
                i++;
            }
        }
        j = strlen;
        if (striptype != LEFTSTRIP)
        {
            j--;
            while (j >= i&&memchr(sep, str[j], charslen))
            {
                j--;
            }
            j++;
        }
        //如果无需要删除的字符
        if (0 == i && j == strlen)
        {
            return str;
        }
        else
        {
            return str.substr(i, j - i);
        }
    }

}

string lstrip(const string & str, const string & chars)
{
    return do_strip( str, LEFTSTRIP, chars);
}


Integer min(const Integer& a, const Integer& b)
{
    if(a <= b) return a;
    return b;
}

Integer max(const Integer& a, const Integer& b)
{
    if(a <= b) return b;
    return a;
}

size_t getMaxBitLength(const Integer& a)
{
    return a.size_bits();
}

EncRectType encryptRectanle(const PaillierFast& crypto, const RectType& rect)
{
    size_t dim = rect[0].size();
    EncRectType encRect;
    for(auto& p : rect){
        encRect.push_back({});
        for(size_t i = 0; i < dim; ++i){
            encRect.back().push_back(crypto.encrypt(p[i]));
        }
    }
    return encRect;
}

RectType decryptRectangle(const PaillierFast& crypto, const EncRectType& encRect)
{
    size_t dim = encRect[0].size();
    RectType rect;
    for(auto& ep: encRect){
        rect.push_back({});
        for(size_t i = 0; i < dim; ++i){
            rect.back().push_back(crypto.decrypt(ep[i]));
        }
    }
    return rect;
}

RectType getMBR(const vector<RectType>& rects)
{
    size_t dim = rects[0][LEFT_BOTTOM_CORNER].size();
    RectType mbr(2);
    for(auto& p: mbr) p.resize(dim);
    for(size_t i = 0; i < dim; ++i){
        mbr[LEFT_BOTTOM_CORNER][i] = std::numeric_limits<int64_t>::max();
        mbr[RIGHT_UP_CORNER][i] = std::numeric_limits<int64_t>::min();
    }
    for(auto& rect: rects){
        for(size_t i = 0; i < dim; ++i){
            mbr[LEFT_BOTTOM_CORNER][i] = min(mbr[LEFT_BOTTOM_CORNER][i], rect[LEFT_BOTTOM_CORNER][i]);
            mbr[RIGHT_UP_CORNER][i] = max(mbr[RIGHT_UP_CORNER][i], rect[RIGHT_UP_CORNER][i]);
        }
    }
    return mbr;
}


ERTreeEntry::ERTreeEntry(const PaillierFast& crypto, const RTreeEntry& entry)
{
    this->level_ = 0;
    this->count_ = crypto.encrypt(0).data;
    this->score_ = crypto.encrypt(0);
    this->is_leaf_ = entry.is_leaf;
    if(entry.is_leaf == false){
        this->rect_ = encryptRectanle(crypto, entry.rect);
        size_t dim = entry.rect[LEFT_BOTTOM_CORNER].size();
        this->midPoint_.resize(dim);
        for(size_t i = 0; i < dim; ++i){
            this->midPoint_[i] = crypto.encrypt(Integer((entry.rect[LEFT_BOTTOM_CORNER][i] + entry.rect[RIGHT_UP_CORNER][i]) / 2));
        }
    }else{
        this->data_.resize(entry.data.size());
        for(size_t i = 0; i < entry.data.size(); ++i){
            this->data_[i] = crypto.encrypt(entry.data[i]);
        }
    }
}

RTree deSeriRTree(const string& filePath)
{
    std::ifstream is(filePath);
    std::stringstream ss;
    ss << is.rdbuf();
    string content(ss.str());
    vector<string> lines;
    splitString(content, lines, "e");
    size_t dim, height;
    {
        vector<string> infos;
        splitString(lines[0], infos, " ");
        dim = std::stoul(infos[0]);
        height = std::stoul(infos[1]);
    }
    shared_ptr<RTreeEntry> root(new RTreeEntry);
    vector<vector<decltype(root)>> levelStruct;
    for(size_t i = 0; i < lines.size() - 1; ++i){
        string level = lstrip(lines[1 + i]);
        if(level.empty()) break;
        levelStruct.push_back({});

        // leaves
        if(i == height - 1){
            // get points
            {
                vector<string> data;
                splitString(level, data, " ");
                for(size_t idx = 1; idx < data.size(); ++idx){
                    // std::cout << data[idx] << std::endl;
                    if(data[idx] == "\n" || data[idx].empty()) continue;
                    vector<string> vals;
                    splitString(data[idx], vals, ",");
                    shared_ptr<RTreeEntry> entry(new RTreeEntry);
                    entry->is_leaf = true;
                    entry->data.resize(vals.size());
                    for(size_t j = 0; j < vals.size(); ++j){
                        entry->data[j] = std::stoul(vals[j]);
                    }
                    levelStruct.back().push_back(std::move(entry));
                }
            }
        }else{
            vector<string> items;
            splitString(level, items, " ");
            size_t size = std::stoul(items[0]);
            for(size_t j = 0; j < size; ++j){
                shared_ptr<RTreeEntry> entry(new RTreeEntry);
                entry->is_leaf = false;
                // get rect
                {
                    entry->rect.resize(2);
                    for(auto& p: entry->rect) p.resize(dim);
                    vector<string> leftBottomCornerVals;
                    vector<string> rightUpCornerVals;
                    splitString(items[1 + j * 3], leftBottomCornerVals, ",");
                    splitString(items[1 + j * 3 + 1], rightUpCornerVals, ",");
                    // std::cout << items[1 + j * 3] << std::endl;
                    // std::cout << items[1 + j * 3 + 1] << std::endl;

                    for(size_t k = 0; k < dim; ++k){
                        entry->rect[LEFT_BOTTOM_CORNER][k] = std::stoul(leftBottomCornerVals[k]);
                        entry->rect[RIGHT_UP_CORNER][k] = std::stoul(rightUpCornerVals[k]);
                    }
                }
                entry->entries.resize(std::stoul(items[1 + j * 3 + 2]));
                levelStruct.back().push_back(std::move(entry));
            }
        }
        if(levelStruct.size() > 1){
            size_t idx = 0;
            for(size_t j = 0; j < levelStruct[levelStruct.size() - 2].size(); ++j){
                for(size_t k = 0; k < levelStruct[levelStruct.size() - 2][j]->entries.size(); ++k){
                    levelStruct[levelStruct.size() - 2][j]->entries[k] = levelStruct.back()[idx];
                    idx += 1;
                }
                
            }
        }
    }
    RTree t;
    if(levelStruct[0].size() > 1){
        shared_ptr<RTreeEntry> root(new RTreeEntry);
        root->is_leaf = false;
        vector<RectType> rects(levelStruct[0].size());
        for(size_t i = 0; i < rects.size(); ++i){
            rects[i] = levelStruct[0][i]->rect;
        }
        root->rect = getMBR(rects);
        root->entries = levelStruct[0];
        t.root = root;
        t.height = height + 1;
        return t;
    }
    t.root = levelStruct[0][0];
    t.height = height;
    return t;
}

shared_ptr<ERTreeEntry> encryptRTree(const PaillierFast& crypto, RTree& t)
{
    int height = t.getHeight();
    shared_ptr<ERTreeEntry> root(new ERTreeEntry(crypto, *(t.root)));//std::make_shared<ERTreeEntry>(crypto, *(t.root));
    
    entryInfos[root.get()] = format("root:[min:%s, max:%s]", point2str(t.root->rect[LEFT_BOTTOM_CORNER]).c_str(), point2str(t.root->rect[RIGHT_UP_CORNER]).c_str()) + ", score:%s";
    root->level_ = height;

    // level order
    queue<pair<shared_ptr<RTreeEntry>, shared_ptr<ERTreeEntry>>> q;
    for(auto entry: t.root->entries){
        q.push({entry, root});
    }
    int i = height - 1;
    while(!q.empty()){
        int size = q.size();
        while(size--){
            auto [curNode, par] = q.front();
            q.pop();
            shared_ptr<ERTreeEntry> eentry(new ERTreeEntry(crypto, *curNode));
            
            if(curNode->is_leaf == false){
                entryInfos[eentry.get()] = format("nonLeaf:[min:%s, max:%s]", point2str(curNode->rect[0]).c_str(), point2str(curNode->rect[1]).c_str()) + ", score:%s";
            }else{
                entryInfos[eentry.get()] = format("Leaf:[%s]", point2str(curNode->data).c_str()) + ", score:%s";
            }
            
            eentry->level_ = i;
            par->entries_.push_back(eentry);

            if(curNode->is_leaf == false){
                for(auto child: curNode->entries){
                    q.push({child, eentry});
                }
            }
        }
        i -= 1;
    }
    return root;
}

void traverseERTree(shared_ptr<ERTreeEntry> root, std::function<void(shared_ptr<ERTreeEntry>)>& fn)
{
    fn(root);
    if(root->is_leaf_ == true){
        return;
    }else{
        for(auto entry: root->entries_){
            traverseERTree(entry, fn);
        }
    }
}

void countRecur(shared_ptr<ERTreeEntry> root)
{
    if(root->is_leaf_ == true){
        root->count_ = 1;
    }else{
        Integer acc = 0;
        for(auto entry: root->entries_){
            countRecur(entry);
            acc += entry->count_;
        }
        root->count_ = acc;
    }
}

void count(shared_ptr<ERTreeEntry> root, const PaillierFast& crypto)
{
    countRecur(root);
    std::function<void(shared_ptr<ERTreeEntry>)> fn = [&crypto](shared_ptr<ERTreeEntry> entry){
        entry->count_ = crypto.encrypt(entry->count_).data;
    };
    traverseERTree(root, fn);
}

vector<PointType> read(const string& path)
{
    std::ifstream is(path, std::ios::in);
	string line;
	vector<PointType> data;
	while (getline(is, line))
	{
		// 打印整行字符串
		//cout << lineStr << endl;
		if(line.empty()) continue;
        vector<string> vals;
		splitString(line, vals, " ");
        data.emplace_back(vals.size());
		for(size_t j = 0; j < vals.size();++j){
			data.back()[j] = Integer(vals[j].c_str());
		}
	}
	is.close();
	return data;
}

vector<shared_ptr<ERTreeEntry>>& B(shared_ptr<ERTreeEntry> entry)
{
    assert(entry->level_ == 2);
    return entry->entries_;
}