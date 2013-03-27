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

            SimpaticoWorkItem work;
            bool done = false;
            while (! done)
                {
                try
                    {
                    work = queue.wait_and_pop();
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
           
            Scalar3 L = box.getL();
            Util::Vector lengths(L.x, L.y, L.z);
            McMd::MdSystem *system_ptr = &simulation->system();
            system_ptr->boundary().setOrthorhombic(lengths);

            McMd::System::MoleculeIterator molIter;
            McMd::Molecule::AtomIterator atomIter;
            unsigned int tag = 0;
            for (int iSpec = 0; iSpec < simulation->nSpecies(); ++iSpec)
                {
                system_ptr->begin(iSpec, molIter);
                for ( ; !molIter.isEnd(); ++molIter)
                    {
                    for (molIter->begin(atomIter); !atomIter.isEnd(); ++atomIter)
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


Simpatico::Simpatico(boost::shared_ptr<SystemDefinition> sysdef,
                     boost::python::object callback,
                     const unsigned int queue_limit)
    : Analyzer(sysdef),m_callback(callback),
      m_work_queue(queue_limit/sizeof(SimpaticoWorkItem))
    {
    }

Simpatico::~Simpatico()
    {
    }

void Simpatico::resetStats()
    {
    // Get parameter file contents from python callback
    boost::python::object rv = m_callback();
    boost::python::extract<std::string> extracted_rv(rv);
    std::string params;
    if (extracted_rv.check())
        {
        params = extracted_rv();
        }
    else
        {
        cerr << "***Error! Parameter file generation failed." << endl;
        throw runtime_error("Cannot initialize simpatico.");
        }

#ifdef ENABLE_MPI
    if (m_pdata->getDomainDecomposition())
        if (! (m_exec_conf->getRank() == 0))
            return;
#endif

    // create one worker thread
    m_worker_thread = boost::thread(SimpaticoWorkerThread(), boost::ref(m_work_queue), params);
    }

void Simpatico::analyze(unsigned int timestep)
    {
    if (m_prof)
        m_prof->push("Simpatico");

    SnapshotParticleData snap(m_pdata->getNGlobal());
    m_pdata->takeSnapshot(snap);

#ifdef ENABLE_MPI
    if (m_pdata->getDomainDecomposition())
        if (! (m_exec_conf->getRank() == 0))
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
    if (m_pdata->getDomainDecomposition())
        if (! (m_exec_conf->getRank() == 0))
            return;
#endif

    // finish worker thread
    m_worker_thread.interrupt();
    m_worker_thread.join();

    }

void export_Simpatico()
    {
    class_< Simpatico, boost::shared_ptr<Simpatico>, bases<Analyzer>, boost::noncopyable>
    ("Simpatico", init< boost::shared_ptr<SystemDefinition>, boost::python::object, unsigned int >());
    }
