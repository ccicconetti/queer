/*
              __ __ __
             |__|__|  | __
             |  |  |  ||__|
  ___ ___ __ |  |  |  |
 |   |   |  ||  |  |  |    Ubiquitous Internet @ IIT-CNR
 |   |   |  ||  |  |  |    C++ quantum routing libraries and tools
 |_______|__||__|__|__|    https://github.com/ccicconetti/quantum-routing

Licensed under the MIT License <http://opensource.org/licenses/MIT>
Copyright (c) 2022 C. Cicconetti <https://ccicconetti.github.io/>

Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this software and associated  documentation files (the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell
copies  of  the Software,  and  to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR
IMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include "Support/random.h"

#include <string>
#include <vector>

namespace uiiit {
namespace qr {

class MecQkdWorkload final
{
 public:
  struct AppInfo {
    unsigned long theRegion = 0; //!< region identifier
    double        theWeight = 0; //!< weight to determine occurence probability
    double        theLoad   = 0; //!< load requested, in operations/second
    double        theRate   = 0; //!< QKD rate requested, in b/s

    std::string toString() const;
  };

  /**
   * @brief Construct a new MecQkdWorkload object from its app info data.
   *
   * @param aAppInfo The info to be copied into this data structure.
   * @param aRv The r.v. to draw randomly the apps.
   * @param aWeighted If true then sampling is based on the app's weights.
   *
   * @throw std::runtime_error if no apps are passed.
   */
  MecQkdWorkload(const std::vector<AppInfo>& aAppInfo,
                 support::RealRvInterface&   aRv,
                 const bool                  aWeighted);

  /**
   * @brief Make a new MecQkdWorkload object from a csv file.
   *
   * @param aFilename The csv file where to read the data from.
   * @param aRv The r.v. to draw randomly the apps.
   * @param aWeighted If true then sampling is based on the app's weights.
   *
   * @return MecQkdWorkload
   */
  static MecQkdWorkload fromCsvFile(const std::string&        aFilename,
                                    support::RealRvInterface& aRv,
                                    const bool                aWeighted);

  /**
   * @brief Draw randomly an app.
   */
  AppInfo operator()();

  /**
   * @brief Return the regions in the app infos.
   *
   * @return std::vector<unsigned long>
   */
  const std::set<unsigned long>& regions() const {
    return theRegions;
  }

 private:
  const std::vector<AppInfo> theAppInfo;
  support::RealRvInterface&  theRv;
  const bool                 theWeighted;
  std::set<unsigned long>    theRegions; // never changed after construction
  std::vector<double>        theWeights; // same
};

} // namespace qr
} // namespace uiiit
