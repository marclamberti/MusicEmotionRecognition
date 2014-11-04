#include <iostream>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>


std::string	GetFileExtension(std::string file) {
	std::string::size_type idx = file.rfind('.');
	if (idx != std::string::npos)
		return file.substr(idx + 1);;
	return std::string("");
}

int 	ExecCommand(std::string const &cmd) {
	int		status;
	pid_t	pid;
	
	if ((pid = fork()) < 0) {
			std::cerr << "Fork failed" << std::endl;
		} else if (pid == 0) {
			std::system(cmd.c_str());
			exit(1);
		} else {
			if (waitpid(pid, &status, 0) > 0) {
				if (WIFEXITED(status)) 
					return WIFEXITED(status);
			}
	}
	return 1;
}

/**
 * Normalize the dataset using ffmpeg
 * Each song is converted into wave (PCM) making it mono with a bit rate of 128kbps
 * and a sampling frequency of 44.1kHz
 * Trim the song to recover 60 seconds that are more likely to contain the chorus.
 * NB: We have defined these 60 seconds to start from 00:45:00 to 01:45:00.
 * Ugly implementation
 * @param file to process and directory in which we should put the processed sound.
 */
void	ProcessSound(std::string const &file, std::string const &directory) {
	std::string extension = GetFileExtension(file);
	if (extension == "mp3") {
		std::string filename = file.substr(0, file.length() - (extension.length() + 1));
		std::string cmd = "ffmpeg -i \"" + directory + "/" + file 
						+ "\" -ar 44100 -ac 1 -codec:a libmp3lame -b:a 128k \"" 
						+ directory + "/" + filename + ".wav\"";
		ExecCommand(cmd);
		cmd = "ffmpeg -i \"" + directory + "/" + filename + ".wav" 
						+ "\" -ss 00:00:45 -t 00:01:00 -acodec copy \"" + directory 
						+ "/" + filename + "_60_seconds.wav\"";
		ExecCommand(cmd);
	}
}


bool	normalize_dataset(std::string const &directory) {
	DIR		*dir;
	dirent	*pdir;

	if ((dir = opendir(directory.c_str())) == NULL){
		std::cerr << "Open directory " << directory << " failed" << std::endl;
		return false;
	}

	while ((pdir = readdir(dir))) {
		if (pdir->d_type == DT_REG) {
			 ProcessSound(pdir->d_name, directory);
		}
	}
	closedir(dir);
	return true;
}

int main(int ac, char **av){
	if (ac < 2) {
		std::cout << "Usage:" << av[0] << " <training_directory>" << std::endl;
		return 1;
	}
	normalize_dataset(av[1]);
    return 0;
}