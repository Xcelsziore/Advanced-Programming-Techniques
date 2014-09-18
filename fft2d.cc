// Distributed two-dimensional Discrete FFT transform
// Aiswaria Nair
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include "Complex.h"
#include "InputImage.h"

using namespace std;



void Transform1D(Complex* h, int w, Complex* H, bool flag)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
   // void Inverse1D(Complex* h, int w, Complex* H);
	Complex partialSum ;
    Complex matrixWeight;
    //int width=w;
    //partialSum = Complex(0,0);
    // Calculation of FFT
    if(flag)
    {
    for(int i = 0; i < w; ++i)
    {
          partialSum = Complex(0,0);
        for(int k = 0; k < w; ++k)
        {   
          
           matrixWeight = Complex(cos(2 * M_PI * i * k/w), -sin(2 * M_PI * i * k/w));
           partialSum  = partialSum + matrixWeight * h[k]; 
		   
	    }
    		H[i] = partialSum; 	   
	
    }
    }
    else
    
    { 
        for(int i = 0; i < w; ++i)
        {
             partialSum = Complex(0,0);
		     for(int k = 0; k < w; ++k)
             {   
		        matrixWeight = Complex(cos(2* M_PI * i* k/w), sin(2* M_PI* i * k/w));
	            partialSum = partialSum + matrixWeight * h[k]; 
          
             }
               // partialSum.real= partialSum.real/w;
               // partialSum.imag= partialSum.imag/w;  
	            H[i].real=partialSum.real/w;
    	        H[i].imag=partialSum.imag/w;
        		//H[i] = partialSum; 
        
    	
    }
}
}



/*void Inverse1D(Complex* h, int w, Complex* H)
{

    Complex partialSum ;
  //  Complex matrixWeight;
  
  for(int n = 0; n < w; n++)
  {    
    //partialSum = Complex(0,0);
	for(int k = 0; k < w; k++)
	{
	  partialSum = partialSum + (Complex((cos(2*M_PI*n*k/w)), (sin(2*M_PI*n*k/w))) * h[k]);		
	}
	partialSum.real = partialSum.real;
	partialSum.imag = partialSum.imag;
	
	H[n]=partialSum;
  }
  
}
*/

void Transform2D(const char* inputFN, bool flag) 
{ 
// Do the 2D transform here.
int counter;
//Number of sub tasks within the process  
int numberOfTasks;
//Rank of every CPU
int rank;

//Image details
int imageHeight;
int imageWidth; 

//InputImage image(outputFN);

InputImage image(inputFN);  // Create the helper object for reading the image
Complex *imageDetails; 
//Pointer to original function after completion of 1D transform
Complex *H;     
            

//MPI SetUp
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfTasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Status message_status;  
  
  //Will be using an array which will capture both, sent and received requests separately             
  MPI_Request *ReceiveReqs, *SentReqs;     

  // Step (1) in the comments is the line above.
  // Your code here, steps 2-9
  
  // Starting row
  int rowStart; 
  int rowEnd;                         
  //bool local=flag;
  
 
  //Message Send receive implementation details
   	int sendingCPU;                              
    long msgLength; 
  
  //Data buffer to send and receive messages 
	double *ReceiveBuffer;                       
    double *SendingBuffer;                     
    int startIndex, endIndex;             

    // get Image dimensions and data
    imageHeight = image.GetHeight();          
    imageWidth = image.GetWidth();          
    imageDetails = image.GetImageData();
  
   //Taking care of H memory allocation
    H = new Complex[(imageHeight/numberOfTasks) * imageWidth]; 
 
  // One Dimensonial transform comes here
     rowStart = (imageHeight/numberOfTasks)*rank;
     
	 counter = imageHeight/numberOfTasks;
    
    for (int i = 0; i < counter; ++i)
    {
        Transform1D (imageDetails + (imageWidth * (i+rowStart)),imageWidth, H + imageWidth*i, flag);
    }
  // In principle with the idea of distributing the task
  
 
    // Split up complex and real here
    
    endIndex = 0;
    
	// One message is composed of Ht*Wdth/Tasks and one row will have data buffer for one task
	
   	 msgLength = 2 *(imageWidth/numberOfTasks)* (imageHeight/numberOfTasks);
	
     ReceiveBuffer = new double[2*imageHeight*imageWidth/numberOfTasks];
     
	 SendingBuffer = new double[2*imageHeight*imageWidth/numberOfTasks];

    // Calculating transpose here
    
    for (int column = 0; column < imageWidth; ++column)
    {
        for (int row = 0; row < imageHeight / numberOfTasks; ++row)
        {
            
			{
			SendingBuffer[endIndex]= H[(row * imageWidth) + column].real;
            ++endIndex;
            SendingBuffer[endIndex]= H[(row * imageWidth) + column].imag;
            ++endIndex;
        }
       /* else
        {
        	SendingBuffer[endIndex]= H[(row * imageWidth) + column].real;
            ++endIndex;
            SendingBuffer[endIndex]= H[(row * imageWidth) + column].imag;
            ++endIndex;	
        	
        }*/
            
        }
    }
    
    // Create receive buffer and requests
    SentReqs = new MPI_Request[numberOfTasks];
    ReceiveReqs = new MPI_Request[numberOfTasks - 1];
   
int returnValue;

// RECEIVE DATA-------------------------------------------------------------------
    for(int position = 0; position < (numberOfTasks - 1); ++position)
    {
        returnValue = MPI_Irecv(ReceiveBuffer + position * msgLength, msgLength, MPI_DOUBLE, MPI_ANY_SOURCE,0,MPI_COMM_WORLD,ReceiveReqs+position);
    
	if (returnValue != MPI_SUCCESS)
			{
				cout << "Rank " << rank
						<< " send failed, returnValue " << returnValue << endl;
				MPI_Finalize();
				exit(1);
			}
			
			}


//SENDING DATA---------------------------------------------------------------------
    
	for(int position = 0; position < numberOfTasks; ++position)
    {
        if (position != rank) 
		{
            returnValue = MPI_Isend(SendingBuffer + position*msgLength, msgLength, MPI_DOUBLE,position, 0,MPI_COMM_WORLD,SentReqs+position);   
            
			if (returnValue != MPI_SUCCESS)
			{
				cout << "Rank " << rank
						<< " send failed, returnValue " << returnValue << endl;
				MPI_Finalize();
				exit(1);
			}
			
			 }
    }

    
    startIndex = msgLength * rank;
    
    for(int row = rank *imageHeight/numberOfTasks; row <(rank + 1)* imageHeight/numberOfTasks; ++row)
    {
        for(int column = rank * imageWidth / numberOfTasks; column < (rank + 1)*imageWidth/numberOfTasks; ++column)
        {
            imageDetails[(row * imageWidth) + column] = Complex(SendingBuffer[startIndex], SendingBuffer[startIndex + 1]);
            startIndex += 2;
        }
    }
  
    // Read function starts here
    
    for(int i = 0; i < numberOfTasks - 1; ++i)
    {
        MPI_Wait(ReceiveReqs+i, &message_status);
        startIndex = i*msgLength;
        sendingCPU = message_status.MPI_SOURCE;
        
        
        
        for(int row = rank * imageHeight/numberOfTasks; row < (rank + 1) * imageHeight / numberOfTasks; ++row )
        {
            for(int column = sendingCPU * imageWidth / numberOfTasks; column < (sendingCPU + 1) * imageWidth / numberOfTasks; ++column)
            {
                imageDetails[(row * imageWidth) + column] = Complex(ReceiveBuffer[startIndex],ReceiveBuffer[startIndex + 1]);
                startIndex += 2;
            }
        }
    }

 //Transform begins here
 
    rowStart = (imageHeight / numberOfTasks) * rank;
    
    for (int j = 0; j < imageHeight / numberOfTasks; ++j)
    {
        Transform1D (imageDetails + (imageWidth*(rowStart+j)), imageWidth, H+ (j*imageWidth), flag);
    }
  
  //Memory Management
  
    delete[] ReceiveBuffer; 
	delete[] SendingBuffer;
    delete[] SentReqs; 
	delete[] ReceiveReqs;
  
//Zeroth Node operations begin here--------------------------------

    msgLength = 2 *imageHeight* imageWidth /numberOfTasks;
//RANK ZERO OPERATIONS COME HERE!-------------------------------------------------------
    
    if (rank == 0)
    {
    	 startIndex = 0;
        // Create receive buffer and Requests
        ReceiveBuffer = new double[msgLength * (numberOfTasks)];
        ReceiveReqs = new MPI_Request[numberOfTasks - 1];

 
        for(int i = 0; i < (numberOfTasks - 1); ++i)
        {
            returnValue = MPI_Irecv(ReceiveBuffer+msgLength*i,msgLength, MPI_DOUBLE, MPI_ANY_SOURCE,0, MPI_COMM_WORLD,i+ReceiveReqs);
           
		   if (returnValue != MPI_SUCCESS)
			{
				cout << "Rank " << rank
						<< " send failed, returnValue " << returnValue << endl;
				MPI_Finalize();
				exit(1);
			}
			 
        }
        for(int column = 0; column <imageWidth / numberOfTasks; ++column)
        {
            for(int row = 0; row < imageHeight; ++row)
            {
                imageDetails[imageWidth*row + column] = H[startIndex];
                ++startIndex;
            }
        }
//Work on deciphering the data that has come from all nodes except zeroth node

        for(int index = 0; index < numberOfTasks-1; ++index)
        {
            MPI_Wait(ReceiveReqs + index, &message_status);
            sendingCPU = message_status.MPI_SOURCE;
            startIndex = msgLength * index;
            
            for(int column = sendingCPU * imageWidth/numberOfTasks; column <( 1+ sendingCPU) *imageWidth /numberOfTasks; ++column)
            {
                   for(int row =0; row < imageHeight; ++row)
                {
                    imageDetails[row * imageWidth+column] = Complex(ReceiveBuffer[startIndex], ReceiveBuffer[startIndex +1]);
                    startIndex += 2;
                }
            }
        }
  
       //Memory Management for receiver units
        delete[] ReceiveReqs;
        delete[] ReceiveBuffer;
    }

// FOR ALL OTHER NODES----------------------------------------------------------------------------------
    else
    {
        endIndex = 0;
		SendingBuffer = new double[msgLength];

 //Buffer it first and then send the  data across
        for(int i = 0; i < imageHeight*imageWidth/numberOfTasks; ++i)
        {
			SendingBuffer[endIndex] = H[i].real;
            ++endIndex;
            SendingBuffer[endIndex] = H[i].imag;
            ++endIndex;
        }
       /* else
        
        {
        	SendingBuffer[endIndex] = H[i].real/imageWidth;
            ++endIndex;
            SendingBuffer[endIndex] = H[i].imag/imageWidth;
            ++endIndex;
        	
        	
        }*/
        

        SentReqs = new MPI_Request[0];
        returnValue = MPI_Isend(SendingBuffer, msgLength, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, SentReqs);
        MPI_Wait(SentReqs, &message_status);

//MEMORY MANAGEMENT---------------------------------------------------
        
		delete[] SendingBuffer;
        delete[] SentReqs;
    }  
    delete[] H ; 
  //For Node One, final task of saving the output comes here-------------------------
  //Writing to Output File    
    if(rank == 0)
	{
        
		if(flag)
		{
			image.SaveImageData("TowerOutputDFT.txt", imageDetails, imageWidth, imageHeight);
		}
		else
		{
            image.SaveImageData("TowerOutputIDFT.txt", imageDetails, imageWidth, imageHeight);
		}
    }
   // MPI_Finalize();
 }
// 2DTransform ends here



  
int main(int argc, char** argv)
{
	int returnValue = MPI_Init(&argc,&argv);
	if (returnValue != MPI_SUCCESS) 
	{
		printf ("Error initiating the MPI pgm. Mission aborted..\n");
		MPI_Abort(MPI_COMM_WORLD, returnValue);
	}

	string fn("Tower.txt"); // default file name
 // MPI_Init(&argc, &argv);   
	if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  //2D DFT
	Transform2D(fn.c_str(),true);  // Perform the transform
	
  //IDFT
  // string fn2("TowerOutputDFT.txt"); // another file name
  // if (argc > 1) fn2 = string(argv[1]);
   //Transform2D(fn2.c_str(),false); // Perform the transform   
	 
  MPI_Finalize(); 
    return 0;
}  


  
