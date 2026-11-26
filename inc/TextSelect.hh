#ifndef _TextSelect_hh
#define _TextSelect_hh

#include <ROOT/REveElement.hxx>
#include <mutex>
#include <tuple>

namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
class TextSelect : public ROOT::Experimental::REveElement
{
  public:
    // Constructor remains implicit or default
    TextSelect() : REveElement{"TextSelect"} {} 
      
    // Setter (Called by the GUI thread/EventDisplayManager)
    void set(int run, int subrun, int event);

    // Getter (Called by the art module thread)
    std::tuple<int, int, int> getRunEvent();

    private:
     int runN = 0;
     int subrunN = 0;
     int eventN = 0;
     std::mutex _mutex; // <<< The synchronization tool
};
}

#endif
