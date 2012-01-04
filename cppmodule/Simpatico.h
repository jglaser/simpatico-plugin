#include <hoomd/Analyzer.h>

#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/diagnostics/DiagnosticManager.h>

#include <boost/python.hpp>

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

    DiagnosticSimulation *m_simulation; //!< ptr to Simpatico simulation object
    boost::python::object m_callback; //!< python callback to generate parameter file
    };

void export_Simpatico();
