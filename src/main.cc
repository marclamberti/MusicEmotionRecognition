#include <iostream>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

/**
 * Normalize the dataset using ffmpeg
 * Each song is converted into wave (PCM) making it mono with a bit rate of 128kbps
 * and a sampling frequency of 44.1kHz
 * Ugly implementation
 * @param directory of the training set
 */
bool	normalize_dataset(std::string const &directory) {
	DIR		*dir;
	dirent	*pdir;
	int		status;
	int		value;
	pid_t	pid;

	if ((dir = opendir(directory.c_str())) == NULL){
		std::cerr << "Open directory " << directory << " failed" << std::endl;
		return false;
	}

	while (pdir = readdir(dir)) {
		std::string	file = pdir->d_name;
		std::string::size_type idx;
		if (pdir->d_type == DT_REG) {
			idx = file.rfind('.');
			if (idx != std::string::npos) {
				std::string extension = file.substr(idx + 1);
				if (extension == "mp3") {
					std::string filename = file.substr(0, file.length() - (extension.length() + 1));
					std::string cmd = "ffmpeg -i \"" + directory + "/" + file + "\" -ar 44100 -ac 1 -codec:a libmp3lame -b:a 128k \"" + directory + "/" + filename + ".wav\"";
					if ((pid = fork()) < 0) {
						std::cerr << "Fork failed" << std::endl;
					} else if (pid == 0) {
						std::system(cmd.c_str());
						exit(1);
					} else {
						if (waitpid(pid, &status, 0) > 0) {
							if (WIFEXITED(status)) {
								value = WIFEXITED(status);
							}
						}
					}
				}
			}
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