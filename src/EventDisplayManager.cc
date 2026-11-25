#include "EventDisplay/inc/EventDisplayManager.hh"
namespace mu2e {

// --- Constructor Implementation ---
/**
 * @brief Constructs the EventDisplayManager singleton instance.
 * * It stores pointers to the global REve manager and the synchronization primitives (mutex/cv).
 * The TextSelect pointer is now passed via the setTextSelectId method instead of the constructor.
 * @param eveMgr The global REveManager instance.
 * @param cv The condition variable for thread signaling.
 * @param m The mutex for thread synchronization.
 * @param fGui The pointer to the custom GUI element.
 */
EventDisplayManager::EventDisplayManager(
    ROOT::Experimental::REveManager* eveMgr,
    std::condition_variable& cv,
    std::mutex& m,
    GUI *fGui)
    : 
    REveElement{"EventManager"}, 
    eveMng_{eveMgr}, 
    cv_{&cv}, 
    m_{&m}, 
    fGui_(fGui)
{
}

// --- Thread Synchronization & Commands ---

/*Signals the waiting Art thread to load the next event.
 This function is invoked by the "NextEvent" REve command button.
*/
void EventDisplayManager::NextEvent()
{   std::unique_lock lock{*m_};
    advance_requested_ = true;
    cv_->notify_all();
}

  void EventDisplayManager::EventUpdateDone(){
    {
      std::lock_guard<std::mutex> lock(display_mutex_);
      display_update_done_ = true;
    }
    display_cv_.notify_one();
  }

/*Terminates the application and exits the process.
 Invoked by the "QuitRoot" REve command button.
 */
void EventDisplayManager::QuitRoot()
{
    std::cout<<"Exit Signal 15, leaving REve Display "<<std::endl;
    exit(15);
}



// --- Element ID Management ---

/**
 * @brief Stores the unique REve Element ID (EId_t) for the TextSelect object.
 * * This ID is used for robust, thread-safe lookup of the object via gEve.
 * @param textId The assigned REve Element ID (e.g., 4285).
 */
void EventDisplayManager::setTextSelectId(std::uint32_t textId) {
    fTextId_ = textId;
    std::cout << "[EventDisplayManager::setTextSelectId] fTextId_ set to: " << fTextId_ << std::endl;
}

// --- Core Lookup and Command Execution ---

/**
 * @brief Executes the command to set the user-requested Run/Event ID.
 * * This runs in the REve thread and communicates the command to the TextSelect object.
 */
void mu2e::EventDisplayManager::goToRunEvent(int runId, int subrunId, int eventId)
{
    std::cout << "[EventDisplayManager::goToRunEvent] received: "
              << runId << "/" << subrunId << "/" << eventId << std::endl;

    TextSelect* fText_obj = nullptr;

    if (ROOT::Experimental::gEve != nullptr) {
        ROOT::Experimental::REveElement* element = ROOT::Experimental::gEve->FindElementById(fTextId_);
        if (element != nullptr)
            fText_obj = dynamic_cast<TextSelect*>(element);
    }

    if (fText_obj == nullptr) {
        std::cerr << "CRITICAL ERROR: TextSelect object not found via gEve->FindElementById("
                  << fTextId_ << "). Cannot set Run/Subrun/Event." << std::endl;
        return;
    }

    fText_obj->set(runId, subrunId, eventId);
}

}
