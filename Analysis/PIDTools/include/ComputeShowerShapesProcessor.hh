#ifndef ComputeShowerShapesProcessor_hh
#define ComputeShowerShapesProcessor_hh 1

#include "EVENT/LCCollection.h"
#include <marlin/Processor.h>

#include "ClusterShapes.h"

class ComputeShowerShapesProcessor : public marlin::Processor {
public:
  ComputeShowerShapesProcessor(const ComputeShowerShapesProcessor&) = delete;
  ComputeShowerShapesProcessor& operator=(const ComputeShowerShapesProcessor&) = delete;

  virtual Processor* newProcessor() override { return new ComputeShowerShapesProcessor; }
  ComputeShowerShapesProcessor();
  virtual void init() override;
  virtual void processRunHeader(lcio::LCRunHeader* run) override;
  virtual void processEvent(lcio::LCEvent* evt) override;
  virtual void check(lcio::LCEvent* evt) override;
  virtual void end() override;

private:
  ClusterShapes* pClusterShapes{};
  std::string _PfoCollection{};
  std::string _ClusterCollection{};
  float _X01{}, _X02{};
  float _Rm1{}, _Rm2{};
  lcio::LCCollection* _PFOCol{};
};

#endif
