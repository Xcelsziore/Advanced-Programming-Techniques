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

#include "Complex.h"
#include "InputImage.h"

using namespace std;


void Transform2D(const char* inputFN) 
{ 
// Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  // 4) Obtain a pointer to the Complex 1d array of input data
  // 5) Do the individual 1D transforms on the rows assigned to your CPU
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().

//Number of sub tasks within the process  
int numberOfTasks;
//Rank of every CPU
int rank;

//Image details
int imageHeight;
int imageWidth; 


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
    H = new Complex[(imageHeight / numberOfTasks) * imageWidth]; 
 
  // One Dimensonial transform comes here
    rowStart = (imageHeight / numberOfTasks)*rank;
    int counter= imageHeight/numberOfTasks;
    for (int i = 0; i < counter; ++i)
    {
        Transform1D (imageDetails + (imageWidth * (i+rowStart)),imageWidth, H+imageWidth*i);
    }
  // In principle with the idea of distributing the task
  
  
  
    // Not using an array of Complex, but instead splitting the real and complex parts up (hence the 2 and +=2 that occurs later)
    endIndex = 0;
	// One message is composed of Ht*Wdth/Tasks and one row will have data buffer for one task
	msgLength = 2 * (imageHeight/numberOfTasks)*(imageWidth/numberOfTasks);
     ReceiveBuffer = new double[2*imageHeight*imageWidth/numberOfTasks];
	 SendingBuffer = new double[2*imageHeight*imageWidth/numberOfTasks];

    // Calculating transpose here
    
    for (int column = 0; column < imageWidth; ++column)
    {
        for (int row = 0; row < imageHeight / numberOfTasks; ++row)
        {
            SendingBuffer[endIndex]= H[(row * imageWidth) + column].real;
            ++endIndex;
            SendingBuffer[endIndex]= H[(row * imageWidth) + column].imag;
            ++endIndex;
        }
    }
    
    // Create receive buffer and requests
    SentReqs = new MPI_Request[numberOfTasks];
    ReceiveReqs = new MPI_Request[numberOfTasks - 1];
   
int returnValue;

// RECEIVE DATA-------------------------------------------------------------------
    for(int position = 0; position < (numberOfTasks - 1); ++position)
    {
        returnValue = MPI_Irecv(ReceiveBuffer + position* msgLength, msgLength, MPI_DOUBLE, MPI_ANY_SOURCE,0,MPI_COMM_WORLD,ReceiveReqs+position);
    }


//SENDING DATA---------------------------------------------------------------------
    
	for(int position = 0; position < numberOfTasks; ++position)
    {
        if (position != rank) 
		{
            returnValue = MPI_Isend(SendingBuffer + position*msgLength, msgLength, MPI_DOUBLE,position, 0,MPI_COMM_WORLD,SentReqs+position);   
        }
    }

    // Get the data from me to myself instead of using MPI
    startIndex = msgLength * rank;
    
    for(int row = rank *imageHeight/numberOfTasks; row <(rank + 1)* imageHeight/numberOfTasks; ++row)
    {
        for(int column = rank * imageWidth / numberOfTasks; column < (rank + 1)*imageWidth / numberOfTasks; ++column)
        {
            imageDetails[(row * imageWidth) + column] = Complex(SendingBuffer[startIndex], SendingBuffer[startIndex + 1]);
            startIndex += 2;
        }
    }
  
    // Read all the sent data
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

    // Do the transform:
    rowStart = (imageHeight / numberOfTasks) * rank;
    for (int j = 0; j < imageHeight / numberOfTasks; ++j)
    {
        Transform1D (imageDetails + (imageWidth*(rowStart+j)), imageWidth, H+j*imageWidth);
    }
  
    // free up the buffer and request memory
    delete[] ReceiveBuffer;
    delete[] SendingBuffer;
    delete[] SentReqs;
    delete[] ReceiveReqs;
  
    // Send all the data back to rank zero:
    msgLength = 2 * imageWidth *imageHeight/numberOfTasks;
//RANK ZERO OPERATIONS COME HERE!-------------------------------------------------------
    
    if (rank == 0)
    {
    	 startIndex = 0;
        // Create receive buffer and Requests
        ReceiveBuffer = new double[msgLength * (numberOfTasks)];
        ReceiveReqs = new MPI_Request[numberOfTasks - 1];

 
        for(int i = 0; i < (numberOfTasks - 1); ++i)
        {
            returnValue = MPI_Irecv(ReceiveBuffer + msgLength*i,msgLength,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,i+ReceiveReqs);
        }
        for(int column = 0; column <imageWidth / numberOfTasks; ++column)
        {
            for(int row = 0; row < imageHeight; ++row)
            {
                imageDetails[imageWidth*row + column] = H[startIndex];
                ++startIndex;
            }
        }

        // read the recieve data
        for(int index = 0; index < numberOfTasks-1; ++index)
        {
            MPI_Wait(ReceiveReqs + index, &message_status);
            
            sendingCPU = message_status.MPI_SOURCE;
            
            startIndex = msgLength * index;
            
            for(int column = sendingCPU * imageWidth / numberOfTasks; column < (sendingCPU + 1) * imageWidth / numberOfTasks; ++column)
            {
                   for(int row = 0; row < imageHeight; ++row)
                {
                    imageDetails[row * imageWidth + column] = Complex(ReceiveBuffer[startIndex], ReceiveBuffer[startIndex + 1]);
                    startIndex += 2;
                }
            }
        }
  
        // Delete/free data:
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

        SentReqs = new MPI_Request[0];
        returnValue = MPI_Isend(SendingBuffer, msgLength, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, SentReqs);
        MPI_Wait(SentReqs, &message_status);

//MEMORY MANAGEMENT---------------------------------------------------
        delete[] SendingBuffer;
        delete[] SentReqs;
    }
      delete[] H;
  //For Node One, final task of saving the output comes here-------------------------    
    if(rank == 0){
        image.SaveImageData("AishTowerOutput.txt", imageDetails, imageWidth, imageHeight);
    }
 }
// 2DTranform ends here

void Transform1D(Complex* h, int w, Complex* H)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
    
	Complex partialSum ;
    Complex matrixWeight;
    
    // Calculation of FFT
    for(int i = 0; i < w; ++i)
    {
        partialSum = Complex(0,0);
        
        for(int k = 0; k < w; ++k)
        {
            matrixWeight = Complex(cos(2 * M_PI * i * k/w), -sin(2 * M_PI * i * k/w));
            partialSum + = matrixWeight * h[k];
        }
        H[i] = partialSum;
    }
	
}

  
int main(int argc, char** argv)
{
	int returnValue = MPI_Init(&argc,&argv);
	if (returnValue != MPI_SUCCESS) 
	{
		printf ("Error initiating the MPI pgm. Mission aborted..\n");
		MPI_Abort(MPI_COMM_WORLD, returnValue);
	}

	string fn("Tower.txt"); // default file name
	MPI_Init(&argc, &argv);   
	if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
	Transform2D(fn.c_str()); // Perform the transform.
    MPI_Finalize(); 
    return 0;
}  


  
