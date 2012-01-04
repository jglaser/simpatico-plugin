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
//   ParticleDataArrays arrays = m_pdata->acquireReadWrite();
 //   std::cout << arrays.nparticles << std::endl;
    std::cout << m_pdata->getN() << std::endl;
    const BoxDim& box = m_pdata->getBox();
    std::cout <<& m_pdata->getBox() << " " << m_pdata << std::endl;
   
    std::cout <<  box.xhi << " " <<  box.xlo  << " " << box.yhi << std::endl;
 
  //  m_pdata->release();
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
    std::cout << box.xhi << " " <<  box.xlo  << box.yhi << std::endl;
    Util::Vector lengths(box.xhi-box.xlo,box.yhi-box.ylo,box.zhi-box.zlo);
    McMd::MdSystem system = m_simulation->system();
    system.boundary().setLengths(lengths);
    std::cout << lengths;

    ParticleDataArraysConst arraysConst = m_pdata->acquireReadOnly();
    McMd::System::MoleculeIterator molIter;
    McMd::Molecule::AtomIterator atomIter;
    for (int iSpec = 0; iSpec < m_simulation->nSpecies(); ++iSpec)
        {
        system.begin(iSpec, molIter);
        for ( ; !molIter.atEnd(); ++molIter)
            {
            for (molIter->begin(atomIter); !atomIter.atEnd(); ++atomIter)
                {
                unsigned int idx = arraysConst.rtag[atomIter->id()];
                atomIter->position() = Util::Vector(arraysConst.x[idx]+lengths[0]/2.,
                                              arraysConst.y[idx]+lengths[1]/2.,
                                              arraysConst.z[idx]+lengths[2]/2.);
                }
            }
        } 
    m_pdata->release();

    system.pairPotential().buildPairList();

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
