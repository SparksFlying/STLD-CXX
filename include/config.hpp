#pragma once
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <array>
#include <queue>
#include <map>
#include <set>
#include <bitset>
#include <iostream>
#include <memory>
#include "rest_rpc.hpp"
#include <chrono>
#include <thread>
#include "ophelib/paillier_fast.h"
#include "ophelib/vector.h"
#include "ophelib/packing.h"
#include "ophelib/util.h"
#include "ophelib/ml.h"
#include "ophelib/random.h"
#include <cmath>
#include <sstream>
#include <mutex>
#include <atomic>
#include <chrono>
using std::string;
using std::vector;
using std::pair;
using std::array;
using std::shared_ptr;
using std::queue;
using std::priority_queue;
using std::map;
using std::set;
using std::bitset;
using namespace ophelib;
using namespace rest_rpc;
using namespace std::chrono;

namespace config{
    const string workspaceFolder = "/home/zanglang/source/STLD-CXX";
    const string logFolder = workspaceFolder + "/log";
    const string dataFolder = workspaceFolder + "/data";
    const string keysFolder = workspaceFolder + "/keys";
    constexpr uint32_t KEY_SIZE = 1024;

    constexpr size_t FLOAT_EXP = 30;
    constexpr size_t NUM_THREADS = 6;

    const string DO_IP = "127.0.0.1";
    const string DSP_IP = "127.0.0.1";
    const string DAP_IP = "127.0.0.1";
    const string CA_IP = "127.0.0.1";
    
    const int DO_PORT = 10011;
    const int DSP_PORT = 10012;
    const int DAP_PORT = 10013;
    const int CA_PORT = 10004;
};