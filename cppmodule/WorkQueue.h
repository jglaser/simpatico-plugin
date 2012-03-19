// implemented after
// http://www.justsoftwaresolutions.co.uk/threading/implementing-a-thread-safe-queue-using-condition-variables.html

#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>
#include <queue>

template<typename Data>
class WorkQueue
{
    public:
        //! Add work to the queue
        void push(Data const& data)
        {
            boost::mutex::scoped_lock lock(m_mutex);
            m_queue.push(data);
            lock.unlock();
            m_condition_variable.notify_one();
        }

        bool empty() const
        {
            boost::mutex::scoped_lock lock(m_mutex);
            return m_queue.empty();
        }

        bool try_pop(Data& popped_value)
        {
            boost::mutex::scoped_lock lock(m_mutex);
            if(m_queue.empty())
            {
                return false;
            }
            
            popped_value=m_queue.front();
            m_queue.pop();
            return true;
        }

        void wait_and_pop(Data &popped_value)
        {
            boost::mutex::scoped_lock lock(m_mutex);
            while(m_queue.empty())
            {
                m_condition_variable.wait(lock);
            }
            
            popped_value=m_queue.front();
            m_queue.pop();
        }

    private:
        std::queue<Data> m_queue;
        mutable boost::mutex m_mutex;
        boost::condition_variable m_condition_variable;
};
