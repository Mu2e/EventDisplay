#ifndef _EventDisplayManager_hh
#define _EventDisplayManager_hh

#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <iostream>
#include <map>
#include <string>
#include "EventDisplay/inc/GUI.hh"
#include "EventDisplay/inc/TextSelect.hh"
#include "nlohmann/json.hpp"
#include <ROOT/REveManager.hxx>
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

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
    private:
        // Default constructor is private to prevent ROOT from creating uninitialized instances.
        EventDisplayManager() = default;

    public:
        /**
         * @brief Primary constructor for initializing thread synchronization and manager pointers.
         * @param eveMgr The global REveManager instance.
         * @param cv The condition variable for inter-thread signaling (notify/wait).
         * @param m The mutex for protecting access and for condition variable use.
         * @param fGui Pointer to the custom GUI element instance.
         * @param geometry The initial geometry type (e.g., "2020", "2025").
         */
        explicit EventDisplayManager(ROOT::Experimental::REveManager* eveMgr,
                                    std::condition_variable& cv,
                                    std::mutex& m,
                                    GUI *fGui,
                                    const std::string& geometry = "2025");

        // --- REve Command Methods (Invoked by browser buttons) ---

        /**
         * @brief Command to signal the Art analysis thread to load the next event.
         */
        void NextEvent(); 
        
        /**
         * @brief Command to terminate the ROOT application and the job gracefully.
         */
        void QuitRoot();

        void autoplay(int x);

        /**
         * @brief Command to handle user input and trigger loading a specific Run and Event ID.
         * * This function runs in the REve thread and uses fTextId_ to look up the TextSelect object.
         */
        void goToRunEvent(int runId, int eventId);

        // --- Public Members ---

        // Stores the current Run ID, often used by commands or display logic.
        int run{0}; 

        // Stores the unique REve Element ID (EId_t) of the TextSelect object.
        // This is the robust mechanism for looking up the TextSelect element.
        std::uint32_t fTextId_{0}; 
        
        /**
         * @brief Setter to store the unique REve Element ID after the object is added to the World.
         * @param textId The assigned REve Element ID.
         */
        void setTextSelectId(std::uint32_t textId);
        
        /**
         * @brief Set the current geometry type and update the corresponding TextSelect element ID.
         * @param geomType The geometry type (e.g., "2020", "2025").
         */
        void setGeometry(const std::string& geomType);
        
        /**
         * @brief Get the TextSelect element ID for the current geometry.
         * @return The element ID corresponding to the current geometry type.
         */
        std::uint32_t getTextSelectIdForGeometry() const;
        
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
        
        // Current geometry type (e.g., "2020", "2025")
        std::string currentGeometry_{"2025"};
        
        // Mapping of geometry types to TextSelect element IDs
        static const std::map<std::string, std::uint32_t> geometryToElementIdMap;
        
    };
}

#endif /* EventDisplayManager_h */
