#include "ErrorMessage.h"
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS

clock_t InfoMessage::start;
vector<clock_t> InfoMessage::ends;
string InfoMessage::InputName;
string InfoMessage::error_info = "Correct Execution!! No Errors!!";
vector<string> InfoMessage::warning_info;

void InfoMessage::ErrorMessage(string message, int mode)
{
	if (mode == 0)
	{
		//message = "警告：" + message;
		warning_info.push_back(message);
		return;
	}
	if (mode == 1)
	{
		ends.push_back(clock());
		message = "" + message;
		error_info = message;
		EndInfo();
		exit(0);
	}
}

void InfoMessage::EndInfo()
{
	char buffer_d[256];
	char buffer_t[256];

	///////////////////// 获取程序执行日期 ////////////////////////
#ifdef _WIN32||defined _WIN64
	time_t local_time = time(nullptr);
#else
	time_t local_time = time(NULL);
#endif
	tm time_tm = *localtime(&local_time);
	strftime(buffer_t, sizeof(buffer_t), "[ %Y-%m-%d %X ]", &time_tm);
	//////////////////////////////////////////////////////////////


	///////////////////// 获取工作目录 ///////////////////////////
#ifdef _WIN32||defined _WIN64
	_getcwd(buffer_d, 256);
#else
	getcwd(buffer_d, 256);
#endif
	//////////////////////////////////////////////////////////////

	ofstream LogFile;
	LogFile.open(InputName);
	int size1 = warning_info.size();
	//int size2 = error_info.size();
	LogFile << " ----------------------------------------------------------" << '\n';
	LogFile << "|   MODEC：A MOLTEN-SALT-REACTOR SPECIFIC DEPLETION CODE   |" << '\n';
	LogFile << " ----------------------------------------------------------" << '\n';
	LogFile << '\n';
	LogFile << '\n';
	LogFile << "[ " << buffer_d << " ]" << '\n';
	LogFile << buffer_t << '\n';
	LogFile << '\n';
	LogFile << "Time cost at each period: " << '\n';
	LogFile << "-------------------------------------------------" << '\n';
	LogFile << "   Initialization : " << double(ends[0] - start) / CLOCKS_PER_SEC << "s" << '\n'; // clock()得到的返回值必须除以CLOCKS_PER_SEC之后才能得到运行时间
	LogFile << "   Calculation    : " << double(ends[1] - ends[0]) / CLOCKS_PER_SEC << "s" << '\n';
	LogFile << "   Saving         : " << double(ends[2] - ends[1]) / CLOCKS_PER_SEC << "s" << '\n';
	LogFile << "------------------------------" << '\n';
	LogFile << "   Total          : " << double(ends[3] - start) / CLOCKS_PER_SEC << "s" << '\n';
	LogFile << "-------------------------------------------------" << '\n';
	LogFile << '\n';
	LogFile << '\n';
	LogFile << '\n';
	LogFile << "Warning Messages during excuting MODEC: " << '\n';
	LogFile << "-------------------------------------------------" << '\n';
	if (size1 > 0)
	{
		for (int i = 0; i < size1; ++i)
		{
			LogFile << "[" << i+1 << "]: ";
			LogFile << warning_info[i] <<'\n';
			LogFile << '\n';
		}
	}
	else
	{
		LogFile << "Perfect Execution!! No Warnings!!" << '\n';
	}
	LogFile << "-------------------------------------------------" << '\n';
	LogFile << '\n';
	LogFile << '\n';
	LogFile << '\n';
	LogFile << "Error Messages during excuting MODEC: " << '\n';
	LogFile << "-------------------------------------------------" << '\n';
	LogFile << error_info <<'\n';
	LogFile << "-------------------------------------------------" << '\n';
	LogFile.close();
}
