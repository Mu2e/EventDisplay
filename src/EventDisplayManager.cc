#include "REve/inc/EventDisplayManager.hh"


namespace mu2e {

  EventDisplayManager::EventDisplayManager(
    ROOT::Experimental::REveManager* eveMgr,
    std::condition_variable& cv,
    std::mutex& m,
    REveMu2eGUI *fGui)
    : REveElement{"EventManager"}, eveMng_{eveMgr}, cv_{&cv}, m_{&m}, fGui_(fGui)
  {}

  void
  EventDisplayManager::NextEvent()
  {
    std::cout<<"PROCESSING....PLEASE WAIT"<<std::endl;
    std::unique_lock lock{*m_};
    cv_->notify_all();
  }

  void
  EventDisplayManager::QuitRoot()
  {
    std::cout<<"Exit Signal 15, leaving REve Display "<<std::endl;
    exit(15);
  }
  
   

}
