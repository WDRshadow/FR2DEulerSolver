#include "thread_pool.h"

ThreadPool::ThreadPool(const size_t num_threads) : stop(false), active_tasks(0)
{
    for (size_t i = 0; i < num_threads; ++i)
    {
        workers.emplace_back([this]
        {
            while (true)
            {
                std::function<void()> task;
                {
                    std::unique_lock lock(queue_mutex);
                    condition.wait(lock, [this] { return stop || !tasks.empty(); });

                    if (stop && tasks.empty()) return;

                    task = std::move(tasks.front());
                    tasks.pop();
                    ++active_tasks;
                }
                task();
                {
                    std::lock_guard lock(queue_mutex);
                    --active_tasks;
                    condition.notify_all();
                }
            }
        });
    }
}

void ThreadPool::enqueue(const std::function<void()>& task)
{
    {
        std::lock_guard lock(queue_mutex);
        tasks.push(task);
    }
    condition.notify_one();
}

void ThreadPool::enqueue_bulk(const std::vector<std::function<void()>>& bulk_tasks)
{
    {
        std::lock_guard lock(queue_mutex);
        for (const auto& task : bulk_tasks)
        {
            tasks.push(task);
        }
    }
    condition.notify_all();
}

void ThreadPool::join()
{
    std::unique_lock lock(queue_mutex);
    condition.wait(lock, [this]
    {
        return tasks.empty() && active_tasks.load() == 0;
    });
}

ThreadPool::~ThreadPool()
{
    {
        std::lock_guard lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread& worker : workers)
    {
        if (worker.joinable())
            worker.join();
    }
}
