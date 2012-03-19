#include "Simpatico.h"

#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/species/Species.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <util/space/Vector.h>

#include <sstream>

using namespace boost::python;

class SimpaticoWorkerThread
    {
    public:
        // The thread main routine
        void operator ()(WorkQueue<SimpaticoWorkItem>& queue, std::string params)
            {
            simulation = new DiagnosticSimulation;
            init_simulation(params);

            bool done = false;
            SimpaticoWorkItem work;
            while (! done)
                {
                try
                    {
                    queue.wait_and_pop(work);
                    process_work_item(work);
                    }
                catch(boost::thread_interrupted const)
                    {
                    done = true;
                    }
                }
            // finish up work
            while (queue.try_pop(work))
                process_work_item(work);

            simulation->output();

            delete simulation;
            }

        //! Initialize the simulation
        void init_simulation(std::string &params)
            {
            istringstream iss(params);
            simulation->readParam(iss); 

            for (int iSpec = 0; iSpec < simulation->nSpecies(); ++iSpec) {
                McMd::Species *speciesPtr = &simulation->species(iSpec);
                int speciesCapacity = speciesPtr->capacity();

                // Add molecules to system
                for (int iMol = 0; iMol < speciesCapacity; ++iMol)
                    {
                    McMd::Molecule *molPtr = &(speciesPtr->reservoir().pop());
                    simulation->system().addMolecule(*molPtr);
                    }
                }

            simulation->init();
            }

 
        //! Process a work item
        void process_work_item(SimpaticoWorkItem &work_item)
            {
            unsigned int timestep = work_item.timestep;
            SnapshotParticleData& snap = work_item.snapshot;
            BoxDim &box = work_item.box;
            
            Util::Vector lengths(box.xhi-box.xlo,box.yhi-box.ylo,box.zhi-box.zlo);
            McMd::MdSystem *system_ptr = &simulation->system();
            system_ptr->boundary().setLengths(lengths);

            McMd::System::MoleculeIterator molIter;
            McMd::Molecule::AtomIterator atomIter;
            unsigned int tag = 0;
            for (int iSpec = 0; iSpec < simulation->nSpecies(); ++iSpec)
                {
                system_ptr->begin(iSpec, molIter);
                for ( ; !molIter.atEnd(); ++molIter)
                    {
                    for (molIter->begin(atomIter); !atomIter.atEnd(); ++atomIter)
                        {
                        atomIter->position() = Util::Vector(snap.pos[tag].x+lengths[0]/2.,
                                                      snap.pos[tag].y+lengths[1]/2.,
                                                      snap.pos[tag].z+lengths[2]/2.);
                        atomIter->velocity() = Util::Vector(snap.vel[tag].x,
                                                            snap.vel[tag].y,
                                                            snap.vel[tag].z);
                        tag++;
                        }
                    }
                } 

            system_ptr->pairPotential().buildPairList();

            simulation->sample(timestep);
            }

        private:
            DiagnosticSimulation *simulation;
    };


Simpatico::Simpatico(boost::shared_ptr<SystemDefinition> sysdef, boost::python::object callback)
    : Analyzer(sysdef),m_callback(callback) 
    {
    }

Simpatico::~Simpatico()
    {
    }

void Simpatico::resetStats()
    {
#ifdef ENABLE_MPI
    if (m_comm)
        if (m_comm->getMPICommunicator()->rank() != 0)
            return;
#endif

    // Get parameter file contents from python callback
    boost::python::object rv = m_callback();
    boost::python::extract<std::string> extracted_rv(rv);
    std::string params;
    if (extracted_rv.check())
        {
        params = extracted_rv();
        }
    else {
        cerr << "***Error! Parameter file generation failed." << endl;
        throw runtime_error("Cannot initialize simpatico.");
         }

    // create one worker thread
    m_worker_thread = boost::thread(SimpaticoWorkerThread(), boost::ref(m_work_queue), params);
    }

void Simpatico::analyze(unsigned int timestep)
    {
    if (m_prof)
        m_prof->push("Simpatico");

    SnapshotParticleData snap(m_pdata->getNGlobal());
    m_pdata->takeSnapshot(snap);

    if (m_prof)
        m_prof->pop();

#ifdef ENABLE_MPI
    if (m_comm)
        if (m_comm->getMPICommunicator()->rank() != 0)
            return;
#endif

    // post a work item
    m_work_queue.push(SimpaticoWorkItem(timestep, snap, m_pdata->getGlobalBox()));

    if (m_prof)
        m_prof->pop();
    }

void Simpatico::printStats()
    {
#ifdef ENABLE_MPI
    if (m_comm)
        if (m_comm->getMPICommunicator()->rank() != 0)
            return;
#endif

    // finish worker thread
    m_worker_thread.interrupt();
    m_worker_thread.join();

    }

void export_Simpatico()
    {
    class_< Simpatico, boost::shared_ptr<Simpatico>, bases<Analyzer>, boost::noncopyable>
    ("Simpatico", init< boost::shared_ptr<SystemDefinition>, boost::python::object >());
    }
