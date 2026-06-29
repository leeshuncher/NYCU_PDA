#pragma once

#include <chrono>
#include <vector>

#include "types.hpp"

struct partition_result
{
    std::vector<char> part;
    int cut = 0;
    int minSize = 0;
    int maxSize = 0;
};

partition_result partition(parsed_input &parsed, std::chrono::steady_clock::time_point deadline);
