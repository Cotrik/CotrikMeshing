// Log.h
#ifndef __LOG_H__
#define __LOG_H__

#include <fstream>
#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

class Log
{
private:
	Log();
	Log(const Log&);
	Log& operator=(const Log&);

public:
	static Log* Instance();
	void PrintRayIntersectionInfo(std::ofstream &debug, float v, float w, float u, const int number);

private:
	static void Destroy();

	static Log* /*volatile*/ m_pInst;
#ifdef _WIN32
	static HANDLE hMutex;
#else
public:
	static pthread_mutex_t m_lock;
#endif

};
#endif // __LOG_H__
