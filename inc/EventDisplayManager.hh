#ifndef _EventDisplayManager_hh
#define _EventDisplayManager_hh

#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <iostream>
#include "EventDisplay/inc/GUI.hh"
#include "EventDisplay/inc/TextSelect.hh"
#include "nlohmann/json.hpp"
#include <ROOT/REveManager.hxx>

namespace ROOT::Experimental {
  class REveManager;
}

namespace mu2e {

    /**
     * @brief Manages event flow, commands, and synchronization between the Art analysis thread
     * and the dedicated ROOT/REve display thread.
     * Inherits from ROOT::Experimental::REveElement to receive browser commands.
     */
    class EventDisplayManager : public ROOT::Experimental::REveElement {
      std::condition_variable display_cv_;
      std::mutex display_mutex_;
      bool display_update_done_ = false;
    public:
        // Default constructor required by ROOT's dictionary generation mechanism.
        EventDisplayManager() = default; 

        /**
         * @brief Primary constructor for initializing thread synchronization and manager pointers.
         * @param eveMgr The global REveManager instance.
         * @param cv The condition variable for inter-thread signaling (notify/wait).
         * @param m The mutex for protecting access and for condition variable use.
         * @param fGui Pointer to the custom GUI element instance.
         */
        explicit EventDisplayManager(ROOT::Experimental::REveManager* eveMgr,
                                    std::condition_variable& cv,
                                    std::mutex& m,
                                    GUI *fGui);

        // --- REve Command Methods (Invoked by browser buttons) ---

        /**
         * @brief Command to signal the Art analysis thread to load the next event.
         */
        void NextEvent();

      void EventUpdateDone();
        
        /**
         * @brief Command to terminate the ROOT application and the job gracefully.
         */
        void QuitRoot();

        /**
         * @brief Command to handle user input and trigger loading a specific Run and Event ID.
         * * This function runs in the REve thread and uses fTextId_ to look up the TextSelect object.
         */
        void goToRunEvent(int runId, int subrunId, int eventId);

        // --- Public Members ---

        // Stores the current Run ID, often used by commands or display logic.
        int run{0}; 

        // Stores the unique REve Element ID (EId_t) of the TextSelect object.
        // This is the robust mechanism for looking up the TextSelect element.
        std::uint32_t fTextId_{0};

        // Set to true by NextEvent() under m_ before notifying cv_, so analyze()'s
        // cv_.wait() predicate guards against spurious wakeups.
        // Reset to false by analyze() after it wakes up.
        bool advance_requested_{false};
        
        /**
         * @brief Setter to store the unique REve Element ID after the object is added to the World.
         * @param textId The assigned REve Element ID.
         */
        void setTextSelectId(std::uint32_t textId);
        int id_;
        void setid(int id){ id_ = id;}
    private:

        // Pointer to the global REve manager, controlling all visualization.
        ROOT::Experimental::REveManager* eveMng_{nullptr}; 
        
        // Pointer to the condition variable, used to unblock the Art thread (e.g., in analyze()).
        std::condition_variable* cv_{nullptr};             
        
        // Pointer to the mutex, used for thread synchronization (locking data access and cv_ usage).
        std::mutex* m_{nullptr};

        // Pointer to the custom GUI element instance.
        GUI *fGui_{nullptr};                               
        
        // Raw pointer to the TextSelect element. While fTextId_ is preferred for lookup, 
        // this is kept for direct access if necessary (but is less robust).
        TextSelect *fText_{nullptr};
        
    };
}

#endif /* EventDisplayManager_h */
