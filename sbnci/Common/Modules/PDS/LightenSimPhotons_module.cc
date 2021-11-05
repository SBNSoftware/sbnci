/**
 * @file   sbnci/Modules/PDS/LightenSimPhotons_module.cc
 * @brief  Converts a `sim::SimPhotons` data product into `sim::SimPhotonsLite`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 13, 2021
 */

// framework libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "lardataobj/Simulation/SimPhotons.h" // also sim::SimPhotonsLite

// C++ standard library
#include <vector>
#include <string>


// -----------------------------------------------------------------------------
namespace sbn {   class LightenSimPhotons; }

/**
 * @brief Converts a `sim::SimPhotons` data product into `sim::SimPhotonsLite`.
 * 
 * This module can convert one or more `sim::SimPhotons` data products
 * from one producer into matching `sim::SimPhotonsLite` data products.
 * 
 * The time conversion is hard-coded in the same way as it is in the legacy
 * LArG4 code `larg4::OpFastScintillation`.
 * 
 * 
 * 
 * Configuration
 * --------------
 * 
 * * `InputLabel` (input tag, mandatory): the input tag for the module producing
 *   the data product(s) to be converted.
 * * `InputInstances` (list of strings, optional): if specified, the instance
 *   name in `InputLabel` is ignored and the instance names in this parameter
 *   are used instead.
 * 
 * 
 * Output
 * -------
 * 
 * * a single `std::vector<sim::SimPhotonLite>` collection for each of the input
 *   collections; if multiple instance names are specified (`InputInstances`),
 *   each input is converted in its own data product with instance name matching
 *   the input one.
 * 
 * 
 */
class sbn::LightenSimPhotons: public art::SharedProducer {
    
    public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> InputLabel {
      Name("InputLabel"),
      Comment("tag of the module producing the sim::SimPhotons to convert")
      // mandatory
      };
    
    fhicl::Sequence<std::string> InputInstances {
      Name("InputInstances"),
      Comment("if not empty, use these instance names for input tags"),
      std::vector<std::string>{}
      };
    
  }; // struct Config
  
  using Parameters = art::SharedProducer::Table<Config>;
  
  
  LightenSimPhotons(Parameters const& config, art::ProcessingFrame const&);
  
  
  void produce(art::Event& event, art::ProcessingFrame const&) override;
  
  
    private:
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  /// List of data products to convert.
  std::vector<art::InputTag> const fInputTags;
  
  // --- END ---- Configuration parameters -------------------------------------
  
  
  /// Converts a full collection of `sim::SimPhotons` via `convertPhoton()`.
  std::vector<sim::SimPhotonsLite> convertCollection
    (std::vector<sim::SimPhotons> const& photons) const;
  
  /// Converts a `sim::SimPhotons` object into `sim::SimPhotonLite`.
  sim::SimPhotonsLite convertPhotons(sim::SimPhotons const& photons) const;
  
  /// Converts the time of a photon into the time tick for `sbn::SimPhotonLite`.
  int convertPhotonTime(float time) const;
  
  
  /// Expands a module label and many instance names into a list of tags.
  static std::vector<art::InputTag> createTagList(
    art::InputTag const& moduleLabel, std::vector<std::string> const& instances
    );
  
  
}; // class sbn::LightenSimPhotons


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
sbn::LightenSimPhotons::LightenSimPhotons
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedProducer(config)
  , fInputTags
      { createTagList(config().InputLabel(), config().InputInstances()) }
{
  
  for (art::InputTag const& inputTag: fInputTags) {
    consumes<std::vector<sim::SimPhotons>>(inputTag);
    produces<std::vector<sim::SimPhotonsLite>>(inputTag.instance());
  }
  
  async<art::InEvent>(); // everything is thread-safe, I swear
  
} // sbn::LightenSimPhotons::LightenSimPhotons()


// -----------------------------------------------------------------------------
void sbn::LightenSimPhotons::produce
  (art::Event& event, art::ProcessingFrame const&)
{
  
  for (art::InputTag const& inputTag: fInputTags) {
    
    std::vector<sim::SimPhotonsLite> litePhotons = convertCollection
      (event.getProduct<std::vector<sim::SimPhotons>>(inputTag));
    
    event.put(
      std::make_unique<std::vector<sim::SimPhotonsLite>>(std::move(litePhotons)),
      inputTag.instance()
      );
    
  } // for
  
} // sbn::LightenSimPhotons::produce()


// -----------------------------------------------------------------------------
std::vector<sim::SimPhotonsLite> sbn::LightenSimPhotons::convertCollection
  (std::vector<sim::SimPhotons> const& photons) const
{
  auto lightenPhoton = [this](sim::SimPhotons const& photons)
    { return convertPhotons(photons); };
  
  std::vector<sim::SimPhotonsLite> litePhotons;
  litePhotons.reserve(photons.size());
  std::transform
    (photons.begin(), photons.end(), back_inserter(litePhotons), lightenPhoton);
  return litePhotons;
  
} // sbn::LightenSimPhotons::convertCollection()


// -----------------------------------------------------------------------------
sim::SimPhotonsLite sbn::LightenSimPhotons::convertPhotons
  (sim::SimPhotons const& photons) const
{
  
  sim::SimPhotonsLite lite { photons.OpChannel() };
  for (sim::OnePhoton const& photon: photons)
    ++lite.DetectedPhotons[convertPhotonTime(photon.Time)];
  
  return lite;
} // sbn::LightenSimPhotons::convertPhotons()


// -----------------------------------------------------------------------------
int sbn::LightenSimPhotons::convertPhotonTime(float time) const {
  
  // this conversion is a replica of the one used in `sim::OpFastScintillation`
  // from LArSoft v09_29_00:
  return static_cast<int>(time);
  
} // sbn::LightenSimPhotons::convertPhotonTime()


// -----------------------------------------------------------------------------
std::vector<art::InputTag> sbn::LightenSimPhotons::createTagList(
  art::InputTag const& moduleLabel, std::vector<std::string> const& instances
  )
{
  std::vector<art::InputTag> inputTags;
  
  if (instances.empty()) inputTags.push_back(moduleLabel);
  else {
    for (std::string const& instanceName: instances) {
      inputTags.emplace_back
        (moduleLabel.label(), instanceName, moduleLabel.process());
    } // for
  }
  
  return inputTags;
} // sbn::LightenSimPhotons::createTagList()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::LightenSimPhotons)


// -----------------------------------------------------------------------------
