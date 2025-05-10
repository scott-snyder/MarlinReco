#ifndef ComputeShowerShapesProcessor_hh
#define ComputeShowerShapesProcessor_hh 1

#include <marlin/Processor.h>
#include "EVENT/LCCollection.h"

#include "ClusterShapes.h"

using namespace lcio ;
using namespace marlin ;

class ComputeShowerShapesProcessor : public Processor{
public:

  ComputeShowerShapesProcessor(const ComputeShowerShapesProcessor&) = delete;
  ComputeShowerShapesProcessor& operator=(const ComputeShowerShapesProcessor&) = delete;

  virtual Processor*  newProcessor() override { return new ComputeShowerShapesProcessor ; }
  ComputeShowerShapesProcessor();
  virtual void init(  ) override;
  virtual void processRunHeader( LCRunHeader* run) override;
  virtual void processEvent( LCEvent * evt ) override;
  virtual void check( LCEvent * evt ) override;
  virtual void end() override;
 
private:
  ClusterShapes *pClusterShapes{};
  std::string _PfoCollection{};
  std::string _ClusterCollection{};
  float _X01{},_X02{};
  float _Rm1{},_Rm2{};
  LCCollection* _PFOCol{};
};

#endif
