// -*- C++ -*-
/*! @file
 * @brief Timer support
 *
 * A stopwatch like timer.
 */

#ifndef STOPWATCH_H
#define STOPWATCH_H

#include<sys/time.h>
#include<cstdio>
#include<cstdlib>

using namespace std;
/*! @defgroup timer Timer
 *
 * @ingroup qdp
 *
 * @{
 */
class StopWatch 
{
public:
  //! Constructor
  StopWatch();

  //! Destructor
  ~StopWatch();

  //! Reset the timer
  void reset();

  //! Start the timer
  void start();

  //! Stop the timer
  void stop();

  //! Get time in seconds
  double getTimeInSeconds();

private:
  unsigned long sec;
  unsigned long usec;
  bool startedP;
  bool stoppedP;

  struct timeval t_start;
  struct timeval t_end;
};


StopWatch::StopWatch() 
{
  stoppedP=false;
  startedP=false;
}

StopWatch::~StopWatch() {}

void StopWatch::reset() 
{
  startedP = false;
  stoppedP = false;
}

void StopWatch::start() 
{
  int ret_val;
  ret_val = gettimeofday(&t_start, NULL);
  if( ret_val != 0 ) 
  {
    fprintf(stderr, "Gettimeofday failed in StopWatch::start()" );
    abort();
  }
  startedP = true;
  stoppedP = false;
}

void StopWatch::stop() 
{
  if( !startedP ) 
  { 
    fprintf(stderr, "Attempting to stop a non running stopwatch in StopWatch::stop()" );
    abort();
  }

  int ret_val;
  ret_val = gettimeofday(&t_end, NULL);
  if( ret_val != 0 ) 
  {
    fprintf(stderr, "Gettimeofday failed in StopWatch::end()" );
    abort();
  }
  stoppedP = true;
}
    
double StopWatch::getTimeInSeconds()  
{
  long secs=0;
  long usecs=0;
  if( startedP && stoppedP ) 
  { 
    if( t_end.tv_sec < t_start.tv_sec ) 
    { 
      fprintf(stderr, "Critical timer rollover" );
      abort();
    }
    else 
    { 
      secs = t_end.tv_sec - t_start.tv_sec;

      if( t_end.tv_usec < t_start.tv_usec ) 
      {
	secs -= 1;
	usecs = 1000000;
      }
      usecs += t_end.tv_usec - t_start.tv_usec;
    }
  }
  else 
  {
    fprintf(stderr, "Either stopwatch not started, or not stopped" );
    abort();
  }

  return (double)secs + ((double)usecs / 1e6);
}


 
#endif
