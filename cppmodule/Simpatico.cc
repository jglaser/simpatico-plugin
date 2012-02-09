#include "Simpatico.h"

#include <hoomd/ParticleData.h>

#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/species/Species.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/space/Vector.h>

#include <sstream>

using namespace boost::python;

Simpatico::Simpatico(boost::shared_ptr<SystemDefinition> sysdef, boost::python::object callback)
    : Analyzer(sysdef), m_simulation(0),m_callback(callback) 
    {
    }

Simpatico::~Simpatico()
    {
    if (m_simulation)
        delete m_simulation;
    }

void Simpatico::resetStats()
    {
    if (m_simulation)
        delete m_simulation;
    m_simulation = new DiagnosticSimulation();

    // Get parameter file contents from python callback
    boost::python::object rv = m_callback();
    boost::python::extract<std::string> extracted_rv(rv);
    if (extracted_rv.check())
        {
        istringstream iss(extracted_rv());
        m_simulation->readParam(iss);
        }
    else {
        cerr << "***Error! Parameter file generation failed." << endl;
        throw runtime_error("Cannot initialize simpatico.");
         }

    for (int iSpec = 0; iSpec < m_simulation->nSpecies(); ++iSpec) {
        McMd::Species *speciesPtr = &m_simulation->species(iSpec);
        int speciesCapacity = speciesPtr->capacity();

        // Add molecules to system
        for (int iMol = 0; iMol < speciesCapacity; ++iMol)
            {
            McMd::Molecule *molPtr = &(speciesPtr->reservoir().pop());
            m_simulation->system().addMolecule(*molPtr);
            }
        }

    m_simulation->init();
    }

void Simpatico::analyze(unsigned int timestep)
    {
    // retrieve box information
    const BoxDim& box = m_pdata->getBox();

    Util::Vector lengths(box.xhi-box.xlo,box.yhi-box.ylo,box.zhi-box.zlo);
    McMd::MdSystem *system_ptr = &m_simulation->system();
    system_ptr->boundary().setLengths(lengths);

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_rtag(m_pdata->getGlobalRTags(), access_location::host, access_mode::read);

    McMd::System::MoleculeIterator molIter;
    McMd::Molecule::AtomIterator atomIter;
    for (int iSpec = 0; iSpec < m_simulation->nSpecies(); ++iSpec)
        {
        system_ptr->begin(iSpec, molIter);
        for ( ; !molIter.atEnd(); ++molIter)
            {
            for (molIter->begin(atomIter); !atomIter.atEnd(); ++atomIter)
                {
                unsigned int idx = h_rtag.data[atomIter->id()];
                atomIter->position() = Util::Vector(h_pos.data[idx].x+lengths[0]/2.,
                                              h_pos.data[idx].y+lengths[1]/2.,
                                              h_pos.data[idx].z+lengths[2]/2.);
                }
            }
        } 

    system_ptr->pairPotential().buildPairList();

    m_simulation->sample(timestep);
    }

void Simpatico::printStats()
    {
    m_simulation->output();
    }

void export_Simpatico()
    {
    class_< Simpatico, boost::shared_ptr<Simpatico>, bases<Analyzer>, boost::noncopyable>
    ("Simpatico", init< boost::shared_ptr<SystemDefinition>, boost::python::object >());
    }
