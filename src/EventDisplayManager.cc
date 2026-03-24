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
 * @param geometry The initial geometry type (e.g., "2020", "2025").
 */
std::string drawfilename("EventDisplay/config/guiutils_current.txt");
SimpleConfig drawconfig(drawfilename);


EventDisplayManager::EventDisplayManager(
    ROOT::Experimental::REveManager* eveMgr,
    std::condition_variable& cv,
    std::mutex& m,
    GUI *fGui,
    const std::string& geometry)
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
{
    std::unique_lock lock{*m_}; // Acquire lock on the mutex
    cv_->notify_all();          // Notify the waiting Art thread (in analyze())
}

/*Terminates the application and exits the process.
 Invoked by the "QuitRoot" REve command button.
 */
void EventDisplayManager::QuitRoot()
{
    std::cout<<"Exit Signal 15, leaving REve Display "<<std::endl;
    exit(15);
}

/*Place holder, it will be used for online running - streaming events
 */
void EventDisplayManager::autoplay(int x)
{
    std::cout << "EventManger autoplay() ....... " << x << std::endl;
    
    TextSelect* fText_obj = nullptr;
    
    // Check if the global REveManager instance is available.
    if (ROOT::Experimental::gEve != nullptr) {
        
        // 1. Retrieve the generic REveElement using the global manager (gEve) and the Element ID.
        // Use geometry-based lookup instead of hardcoded ID.
        //std::uint32_t elementId = getTextSelectIdForGeometry();
        ROOT::Experimental::REveElement* element = ROOT::Experimental::gEve->FindElementById(drawconfig.getInt("GUIID")); 
        
        if (element != nullptr) {
            // 2. Safely cast the generic element to the specific TextSelect type.
            fText_obj = dynamic_cast<TextSelect*>(element);
        }
    }
    
    if (fText_obj == nullptr) {
        // CRITICAL ERROR if the object wasn't found or the cast failed.
        std::cerr << "CRITICAL ERROR: TextSelect object not found via gEve->FindElementWithId(" 
                  << fTextId_ << "). Cannot set Run/Event." << std::endl; 
        return; 
    }

    // 3. Command executed: Set the Run/Event numbers in the TextSelect object.
    fText_obj->setAutoplay(x); 
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

/**
 * @brief Set the current geometry type.
 * @param geomType The geometry type (e.g., "2020", "2025").
 */
void EventDisplayManager::setGeometry(const std::string& geomType) {
    currentGeometry_ = geomType;
    std::cout << "[EventDisplayManager::setGeometry] Geometry set to: '" << currentGeometry_ << "'" << std::endl;
}

// --- Core Lookup and Command Execution ---

/**
 * @brief Executes the command to set the user-requested Run/Event ID.
 * * This runs in the REve thread and communicates the command to the TextSelect object.
 */
void mu2e::EventDisplayManager::goToRunEvent(int runId, int eventId)
{
    std::cout << "[EventDisplayManager::goToRunEvent] received: " << runId<<" "<<eventId << std::endl;
    
    TextSelect* fText_obj = nullptr;
    
    // Check if the global REveManager instance is available.
    if (ROOT::Experimental::gEve != nullptr) {
        
        // 1. Retrieve the generic REveElement using the global manager (gEve) and the Element ID.
        // Use geometry-based lookup instead of hardcoded ID.
        //std::uint32_t elementId = getTextSelectIdForGeometry();
        ROOT::Experimental::REveElement* element = ROOT::Experimental::gEve->FindElementById(drawconfig.getInt("GUIID")); 
        
        if (element != nullptr) {
            // 2. Safely cast the generic element to the specific TextSelect type.
            fText_obj = dynamic_cast<TextSelect*>(element);
        }
    }
    
    if (fText_obj == nullptr) {
        // CRITICAL ERROR if the object wasn't found or the cast failed.
        std::cerr << "CRITICAL ERROR: TextSelect object not found via gEve->FindElementWithId(" 
                  << fTextId_ << "). Cannot set Run/Event." << std::endl; 
        return; 
    }

    // 3. Command executed: Set the Run/Event numbers in the TextSelect object.
    fText_obj->set(runId, eventId); 
}

}
