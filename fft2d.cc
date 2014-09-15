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
  int i,r,c;                              // Loop variables
  int rc;                                 // MPI return value
 
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
  
   // allocate space for H
    H = new Complex[(imageHeight / numberOfTasks) * imageWidth];
 
   // Do the 1D transform
    rowStart = (imageHeight / numberOfTasks) * rank;
    for (i = 0; i < imageHeight / numberOfTasks; ++i)
    {
        Transform1D(imgdata + (imgWidth * (i + start_row)), imgWidth, H + imgWidth*i);
    }
  
    // Create send buffers, need to be 1 row for each task, with (Height * Width)/numberOfTasks elements per message
    // Not using an array of Complex, but instead splitting the real and complex parts up (hence the 2 and +=2 that occurs later)
    SendingBuffer = new double[2 * imageHeight * imageWidth / numberOfTasks];
    msgLength = 2 * (imageHeight / numberOfTasks) * (imageWidth / numberOfTasks);
    endIndex = 0;

    // Fill the buffers and transpose at the same time:
    for (c = 0; c < imageWidth; ++c)
    {
        for (r = 0; r < imgHeight / numtasks; ++r)
        {
            sendBuff[tgt] = H[(r * imgWidth) + c].real;
            ++tgt;
            sendBuff[tgt] = H[(r * imgWidth) + c].imag;
            ++tgt;
        }
    }
    
    // Create receive buffer and requests
    ReceiveBuffer = new double[2 * imageHeight*imageWidth / numberOfTasks];
    ReceiveReqs = new MPI_Request[numberOfTasks - 1];
    SentReqs = new MPI_Request[numberOfTasks];

    // Receive the data:
    for(i = 0; i < (numberOfTasks - 1); ++i)
    {
        rc = MPI_Irecv(
            recvBuff + i*message_length,    // Location of receive buffer
            message_length,                 // Length of Message
            MPI_DOUBLE,                     // Data Type
            MPI_ANY_SOURCE,                 // Receive from any source
            0,                              // Tag of zero
            MPI_COMM_WORLD,                 // Communicate with all members
            Rrequests+i);                   // Requests array
    }

    // Send the necessary data:
    for(i = 0; i < numberOfTasks; ++i)
    {
        if (i != rank) {
            rc = MPI_Isend(
                sendBuff + i*message_length,// Location of send buffer
                message_length,             // Length of Message
                MPI_DOUBLE,                 // Data Type
                i,                          // Destination rank
                0,                          // Tag of zero
                MPI_COMM_WORLD,             // Communicate with all members
                Srequests+i);               // Send array
        }
    }

    // Get the data from me to myself instead of using MPI
    startIndex = msgLength * rank;
    for(r = rank * imageHeight / numberOfTasks; r < (rank + 1) * imageHeight / numberOfTasks; ++r)
    {
        for(c = rank * imgWidth / numtasks; c < (rank + 1)*imgWidth / numtasks; ++c)
        {
            imgdata[(r * imgWidth) + c] = Complex(sendBuff[source], sendBuff[source + 1]);
            source += 2;
        }
    }
  
    // Read all the sent data
    for(i = 0; i < numberOfTasks - 1; ++i)
    {
        MPI_Wait(Rrequests+i, &message_status);
        sender = message_status.MPI_SOURCE;
        source = message_length*i;
        for(r = rank * imgHeight/numtasks; r < (rank + 1) * imgHeight / numtasks; ++r )
        {
            for(c = sender * imgWidth / numtasks; c < (sender + 1) * imgWidth / numtasks; ++c)
            {
                imgdata[(r * imgWidth) + c] = Complex(recvBuff[source], recvBuff[source + 1]);
                source += 2;
            }
        }
    }

    // Do the transform:
    rowStart = (imageHeight / numberOfTasks) * rank;
    for (i = 0; i < imageHeight / numberOfTasks; ++i)
    {
        Transform1D(imgdata + (imgWidth * (i + start_row)), imgWidth, H + imgWidth * i);
    }
  
    // free up the buffer and request memory
    delete[] ReceiveBuffer;
    delete[] SendingBuffer;
    delete[] SentReqs;
    delete[] ReceiveReqs;
  
    // Send all the data back to rank zero:
    msgLength = 2 * imageWidth * imageHeight / numberOfTasks;

    // If I'm rank zero:
    if (rank == 0)
    {
        // Create receive buffer and Requests
        recvBuff = new double[message_length * (numtasks)];
        Rrequests = new MPI_Request[numtasks - 1];

        // Recieve data:
        for(i = 0; i < numtasks - 1; ++i)
        {
            rc = MPI_Irecv(recvBuff + i * message_length, message_length,  MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, Rrequests + i);
        }
        // fill the data I should send to myself rather than use a message
        source = 0;
        for(c = 0; c < imgWidth / numtasks; ++c)
        {
            for(r = 0; r < imgHeight; ++r)
            {
                imgdata[r * imgWidth + c] = H[source];
                ++source;
            }
        }

        // read the recieve data
        for(i = 0; i < numtasks - 1; ++i)
        {
            MPI_Wait(Rrequests + i, &message_status);
            sender = message_status.MPI_SOURCE;
            source = message_length * i;
            for(c = sender * imgWidth / numtasks; c < (sender + 1) * imgWidth / numtasks; ++c)
            {
                for(r = 0; r < imgHeight; ++r)
                {
                    imgdata[r * imgWidth + c] = Complex(recvBuff[source], recvBuff[source + 1]);
                    source += 2;
                }
            }
        }
  
        // Delete/free data:
        delete[] Rrequests;
        delete[] recvBuff;
    }

    // If I'm not node zero:
    else
    {
        // Create send buffer:
        sendBuff = new double[message_length];
        tgt = 0;

        // Fill up the send buffer with data:
        for(i = 0; i < imgWidth * imgHeight / numtasks; ++i)
        {
            sendBuff[tgt] = H[i].real;
            ++tgt;
            sendBuff[tgt] = H[i].imag;
            ++tgt;
        }
        // Send the data:
        Srequests = new MPI_Request[0];
        rc = MPI_Isend(sendBuff, message_length, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, Srequests);
        MPI_Wait(Srequests, &message_status);

        // Delete/free data:
        delete[] sendBuff;
        delete[] Srequests;
    }

    // Delete the H array:
    delete[] H;

    // Save the image data into the file Tower=DFT2D.txt:
    if(rank == 0)
    {
        image.SaveImageData("Tower-DFT2D.txt", imgdata, imgWidth, imgHeight);
    }
} // End Transform2D


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
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
}  
  
int main(int argc, char** argv)
{
	int rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	string fn("Tower.txt"); // default file name
	if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
	Transform2D(fn.c_str()); // Perform the transform.
	return 0;
}  


  
