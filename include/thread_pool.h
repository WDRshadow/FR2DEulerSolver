#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>

class ThreadPool
{
public:
    explicit ThreadPool(size_t num_threads);
    ~ThreadPool();

    // 添加一个任务
    void enqueue(const std::function<void()>& task);

    // 批量添加任务
    void enqueue_bulk(const std::vector<std::function<void()>>& bulk_tasks);

    // 等待所有任务完成
    void join();

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<bool> stop;

    std::atomic<int> active_tasks; // 正在运行的任务数
};


#endif //THREAD_POOL_H
