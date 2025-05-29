#include <gtest/gtest.h>

#include "thread_pool.h"

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(THREAD_POOL, test)
{
    ThreadPool pool(4);

    std::vector<std::function<void()>> tasks;
    tasks.reserve(10);
    for (int i = 0; i < 10; ++i)
    {
        tasks.emplace_back([i]()
        {
            std::cout << "任务 " << i << " 开始\n";
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::cout << "任务 " << i << " 完成\n";
        });
    }
    pool.enqueue_bulk(tasks);
    pool.join();
}
