#include <hoomd/hoomd.h>


#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/diagnostics/DiagnosticManager.h>

#include <stdio.h>
#include <boost/thread.hpp>
#include <boost/python.hpp>

#include "WorkQueue.h"

struct SimpaticoWorkItem
    {
    //! Empty item
    SimpaticoWorkItem()
        : timestep(0),
          box(BoxDim(1.0,1.0,1.0))
        { }

    SimpaticoWorkItem(const unsigned int _timestep,
                      const SnapshotParticleData& _snapshot,
                      const BoxDim& _box)
        : timestep(_timestep),
          snapshot(_snapshot),
          box(_box)
        { }

    unsigned int timestep;
    SnapshotParticleData snapshot;
    BoxDim box;
    };

// simple wrapper around MdSimulation to access diagnosticManager functions

class DiagnosticSimulation : public McMd::MdSimulation
    {
    public:
    DiagnosticSimulation()
        : McMd::MdSimulation()
        {
        }

    // initialize diagnostics
    void init()
        {
        diagnosticManager().initialize();
        }

    // sample one configuration
    void sample(int iStep)
        {
        diagnosticManager().sample(iStep);
        }

    // output statistics
    void output()
        {
        diagnosticManager().output();
        }
    };

class Simpatico : public Analyzer
    {
    public:

        //! Constructor
        //! \param params simpatico parameter file 
        Simpatico(boost::shared_ptr<SystemDefinition>, boost::python::object callback);
        virtual ~Simpatico();

        //! Use simpatico to perform analysis
        virtual void analyze(unsigned int timestep);

        //! Initialize statistcs at beginning of run
        virtual void resetStats();

        //! Output statistics at end of run
        virtual void printStats();

    private:

        boost::python::object m_callback; //!< python callback to generate parameter file
        
        WorkQueue<SimpaticoWorkItem> m_work_queue;

        boost::thread m_worker_thread;    //!< The worker thread
    };

void export_Simpatico();
