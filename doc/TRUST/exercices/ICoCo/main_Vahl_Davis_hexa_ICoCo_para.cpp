#include <Problem.h>
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"

#include "TrioDEC.hxx"
#include <set>
#include <time.h>
#include <ICoCoTrioField.h>
#include <fstream>
#include <string.h>

using namespace MEDCoupling;
using namespace std;
using namespace ICoCo;

//utility methods for synchronizing
//data from the the two trio instance
typedef enum {sync_and,sync_or} synctype;

int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    CommInterface comm;
    set<int> dom_ids;
    dom_ids.insert(0);
    dom_ids.insert(1);
    MPIProcessorGroup dom_group(comm,dom_ids);
    std::string data_file;
    const MPI_Comm* mpicomm=0;
    if (dom_group.containsMyRank()) {
      data_file="PAR_Vahl_Davis_hexa_ICoCo_para"; // name of the data file // + ICoCo couplage
      // Redirection des sorties couplage dans dom.out et dom.err
      if (!freopen("PAR_Vahl_Davis_hexa_ICoCo_para.out","w",stdout)) abort();         // + ICoCo couplage
      if (!freopen("PAR_Vahl_Davis_hexa_ICoCo_para.err","w",stderr)) abort();         // + ICoCo couplage
      mpicomm=dom_group.getComm();
    }
    else
      throw 0;


    Problem* T;
    T = getProblem();
    T->setDataFile(data_file);
    T->setMPIComm((void*)mpicomm);
    T->initialize();

    vector<string> outputnames= T->getOutputFieldsNames();
    for (unsigned int ii=0;ii<outputnames.size();ii++)
      cout<<data_file<<"  Field Output " <<outputnames[ii]<<endl;
    cout<<endl;

    vector<string> inputnames= T->getInputFieldsNames();
    for (unsigned int ii=0;ii<inputnames.size();ii++)
      cout<<data_file<<"  Field Input  " <<inputnames[ii]<<endl;
    cout<<endl;

    bool stop=false; // Does the Problem want to stop ?
    bool ok=true; // Is the time interval successfully solved ?

    TrioField Temperature_field;

    double* save_old_field=0;

    clock_t clock0= clock ();
    int compti=0;

    bool init=true; // first time step ??

    // loop on time steps
    while (!stop)
      {
        compti++;
        clock_t clocki= clock ();
        cout << compti << " CLOCK " << (clocki-clock0)*1.e-6 << endl; 
        ok=false; // Is the time interval successfully solved ?

        // Loop on the time interval tries
        while (!ok && !stop)
          {
            // Compute the first time step length
            double dt=T->computeTimeStep(stop);
            if (stop) // Impossible to solve the next time step, the Problem has given up
              break;
            // Prepare the next time step
            T->initTimeStep(dt);
            {
                  // name in the jdd 
                  T->getInputFieldTemplate("TEMPERATURE_IN_DOM",Temperature_field);
                  // completing the field
                  Temperature_field.set_standalone();
                  int nbcase=Temperature_field.nb_values()*Temperature_field._nb_field_components;
                  for(int i=0;i<nbcase;i++)
                    Temperature_field._field[i]=10;
                  T->setInputField("TEMPERATURE_IN_DOM",Temperature_field);
            }
            clock_t clock_avt_resolution= clock ();
            // Solve the next time step
            ok=T->solveTimeStep();
            clock_t clock_ap_resolution= clock ();
            cout << compti << " TEMPS DE RESOLUTION DU PB (s) :  " << (clock_ap_resolution-clock_avt_resolution)*1.e-6 << endl; 
            init=false;
            if (!ok) // The resolution failed, try with a new time interval.
              T->abortTimeStep();
            else // The resolution was successful, validate and go to the
              // next time step.
              T->validateTimeStep();
          }
        // Stop the resolution if the Problem is stationnary
        bool stat=T->isStationary();
        if (stat)
          stop=true;
      }
    T->terminate();
    delete T;
    if (save_old_field!=0) 
      delete [] save_old_field;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
