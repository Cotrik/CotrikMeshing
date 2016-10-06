#include "Log.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;

Log* Log::m_pInst = NULL;
#ifdef _WIN32
HANDLE Log::hMutex = CreateMutex(NULL, false, NULL);
#else
pthread_mutex_t Log::m_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

const float EPSILON = 0.00001f;

void Log::PrintRayIntersectionInfo(std::ofstream &debug, float v, float w, float u, const int number)
{
#ifdef _WIN32
	WaitForSingleObject(hMutex, INFINITE);
#endif
// 	if (abs(v) < EPSILON) v = 0.0f;
// 	if (abs(w) < EPSILON) w = 0.0f;
// 	if (abs(u) < EPSILON) u = 0.0f;

	debug << "g_i = " << number << ", v = " << v << ", w = " << w << ", u = " << u  << endl;
#ifdef _WIN32
	ReleaseMutex(hMutex);
#endif
}

Log::Log() 
{
	atexit(Destroy);
}

Log* Log::Instance()
{
	if(m_pInst == NULL)
	{
#ifdef _WIN32
		WaitForSingleObject(hMutex, INFINITE);
#else
		pthread_mutex_lock(&m_lock);
#endif
		if(m_pInst == NULL)
		{
#ifndef _WIN32
			if (pthread_mutex_init(&m_lock, NULL) != 0)
			{
			    std::cerr << "\n mutex init failed\n" << std::endl;
			}
#endif
			m_pInst = new Log;
		}
#ifdef _WIN32
		ReleaseMutex(hMutex);
#else
		pthread_mutex_unlock(&m_lock);
#endif
	}
	return m_pInst;
}

void Log::Destroy() 
{
	delete m_pInst;
	m_pInst = NULL;
#ifndef _WIN32
	pthread_mutex_destroy(&m_lock);
#endif
}
