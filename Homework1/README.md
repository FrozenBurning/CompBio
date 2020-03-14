<!--
 * @Description: 
 * @Author: Zhaoxi Chen
 * @Github: https://github.com/FrozenBurning
 * @Date: 2020-03-14 10:55:36
 * @LastEditors: Zhaoxi Chen
 * @LastEditTime: 2020-03-14 11:21:05
 -->
 
# CompBio - Homework1

## Longest Common Subsequence

Developed both on C++ and python, while C++ version for efficiency.

**Author: 陈昭熹**

**Features:**

- 基于动态规划的LCS算法，亦即Needleman Wunsch Algorithm，时间复杂度$O(mn)$
- Python版本是简单的递归实现，且能够输出所有可能答案，用于小规模测试(**题目1**)
- C++版本进行了优化，用于大规模问题的快速求解(**题目2**)

### Shortcut

[Python Script](python/LCS.py)
[C++ Script](cpp/LCS.cpp)
[题目2结果](Answer.txt)

## 1.Prerequisites

### 1.1 OS

Ubuntu 18.04 LTS

### 1.2 Python 3.6

Python3相关版本均可，不支持Python2

## 2. Usage

### 2.1 Python Script

默认求解目标为AAGC和AGT，如下方式运行：
```bash
cd python
python3 LCS.py
```

求解自定义目标seq1和seq2，如下方式进行参数配置：
```bash
python3 LCS.py --seq1 seq1 --seq2 seq2
```

更多参数配置请使用-h命令获取帮助:
```bash
python3 LCS.py -h
```

### 2.2 C++ Script

编译使用的Makefile已经一并提供,若不能直接运行请重新编译。
可以直接运行，查看以LongestCommonSeq.txt作为输入的结果，输出会存放在Answer.txt中

```bash
cd cpp
./main
```