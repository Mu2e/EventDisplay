#ifndef _EventDisplayManager_hh
#define _EventDisplayManager_hh

#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <iostream>
#include "REve/inc/REveMu2eGUI.hh"
#include "nlohmann/json.hpp"

namespace ROOT::Experimental {
  class REveManager;
}

namespace mu2e {

  constexpr auto invalid_event = std::numeric_limits<unsigned>::max();

  class EventDisplayManager : public ROOT::Experimental::REveElement {
    public:
      EventDisplayManager() = default; // ROOT needs a dictionary

      explicit EventDisplayManager(ROOT::Experimental::REveManager* eveMgr,
                                       std::condition_variable& cv,
                                       std::mutex& m,
                                       REveMu2eGUI *fGui);

    void NextEvent();
    void QuitRoot();
    void autoplay_scheduler()
   {
      while (true)
      {
         bool autoplay;
         {
                std::unique_lock<std::mutex> lock{m_mutex};
                if (!m_autoplay)
                {
                // printf("exit thread pre wait\n");
                return;
                }
                if (m_CV.wait_for(lock, m_deltaTime) != std::cv_status::timeout)
                {
                printf("autoplay not timed out \n");
                if (!m_autoplay)
                {
                    printf("exit thread post wait\n");
                    return;
                }
                else
                {
                    continue;
                }
                }
                autoplay = m_autoplay;
         }
         if (autoplay)
         {
                std::cout<<"I'm true!!"<<std::endl;
         }
         else
         {
                return;
         }
      }
   }

   void autoplay(bool x)
   {
      std::cout << "Set autoplay " << x << std::endl;
      static std::mutex autoplay_mutex;
      std::unique_lock<std::mutex> aplock{autoplay_mutex};
      {
         std::unique_lock<std::mutex> lock{m_mutex};

         StampObjProps();
         m_autoplay = x;
         if (m_autoplay)
         {
                if (m_timerThread)
                {
                m_timerThread->join();
                delete m_timerThread;
                m_timerThread = nullptr;
                }
                NextEvent();
                m_timerThread = new std::thread{[this]
                                                { autoplay_scheduler(); }};
         }
         else
         {
                m_CV.notify_all();
         }
      }
   }

   void playdelay(float x)
   {
      printf("playdelay %f\n", x);
      std::unique_lock<std::mutex> lock{m_mutex};
      m_deltaTime = std::chrono::milliseconds(int(x));
      StampObjProps();
      m_CV.notify_all();
   }
    //void InitGuiInfo();

    private:
      ROOT::Experimental::REveManager* eveMng_{nullptr};
      std::condition_variable* cv_{nullptr};
      std::mutex* m_{nullptr};
      bool doneProcessingEvents_{false};
      REveMu2eGUI *fGui_{nullptr};
      
         std::chrono::duration<double> m_deltaTime{1};
   std::thread *m_timerThread{nullptr};
   std::mutex m_mutex;
   std::condition_variable m_CV;
   bool m_autoplay{false};
    };
}

#endif /* EventDisplayManager_h */
