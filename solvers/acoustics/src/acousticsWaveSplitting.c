#include "acoustics.h"

void acousticsWSExchange(acoustics_t *acoustics){
	mesh_t *mesh = acoustics->mesh;
	dlong tag = 101;
	acoustics->o_vtSend.copyTo(acoustics->vtSend);
  MPI_Request *MPIRequests;
  MPIRequests = (MPI_Request *)malloc((mesh->size-1)*2 * sizeof(MPI_Request));
	dlong MPIRequestCounter = 0;
	for(int i = 0; i < mesh->size; i++){
		// Rank i sends, others receive
		if(i == mesh->rank){

			if(acoustics->NComPointsToSendAllRanks){
				for(int j = 0; j < mesh->size; j++){
					dlong sendCount = acoustics->recvCountsArray[mesh->size*i + j];
					if(sendCount){
						
						// TODO: PRECALCULATE ONCE, INSTED OF EVERY TIME.
						dlong sendOffset = 0;
						for(int k = 0; k < j; k++){
							sendOffset += acoustics->recvCountsArray[mesh->size*i + k]*3;
						}
						
						MPI_Isend(acoustics->vtSend+sendOffset, sendCount*3, MPI_DFLOAT, j, tag,
								mesh->comm, &MPIRequests[MPIRequestCounter]);
						MPIRequestCounter++;
					}
				}
			}
		} else{
			dlong recvCount = acoustics->recvCountsArray[mesh->size*i + mesh->rank];
			if(recvCount){
				
				// TODO: PRECALCULATE ONCE, INSTED OF EVERY TIME.
				dlong recvOffset = 0;
				for(int k = 0; k < i; k++){
					recvOffset += acoustics->recvCountsArray[mesh->size*k + mesh->rank]*3;
				}

				MPI_Irecv(acoustics->vtRecv+recvOffset, recvCount*3, MPI_DFLOAT, i, tag, 
						mesh->comm, &MPIRequests[MPIRequestCounter]);
				MPIRequestCounter++;
			}

		}
		
	}	
	MPI_Status *status = (MPI_Status*) calloc(MPIRequestCounter, sizeof(MPI_Status));
  
	MPI_Waitall(MPIRequestCounter, MPIRequests, status);
	free(status);
	acoustics->o_vtRecv.copyFrom(acoustics->vtRecv);
}