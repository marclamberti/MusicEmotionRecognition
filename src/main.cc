#include <iostream>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <aubio/aubio.h>
#include <vector>

namespace {
	std::vector<std::vector<float>> features;

	const int size = 512;
}

bool	getBeat(const std::string &file) {
	aubio_tempo_t 		*tempo;
  	uint_t 				win_size 		= size;
  	uint_t 				hop_size 		= size / 4;
  	uint_t 				n_frames 		= 0;
  	uint_t				read 			= 0;
	aubio_source_t		*aubio_source 	= NULL;
	uint_t				samplerate		= 0;

	// Create a new source object from the file given as parameter
	// If samplerate is set to 0, the sample rate of the original file is used
	// hop_size is the size of the blocks to read from the source
	if (!(aubio_source = new_aubio_source (const_cast<char *>(file.c_str()), samplerate, hop_size))) {
		std::cerr << "new_aubio_source failed" << std::endl;
		aubio_cleanup();
		return false;
	}
	
	// Get the samplerate from the file
	if ((samplerate = aubio_source_get_samplerate(aubio_source)) == 0) {
		std::cerr << "aubio_source_get_samplerate failed" << std::endl;
		del_aubio_source(aubio_source);
		aubio_cleanup();
		return false;
	}
	
	// Create two vectors
	// One input buffer and one for the output position
	fvec_t *in 	= new_fvec(hop_size);
  	fvec_t *out = new_fvec(2);

  	// Create the tempo object
  	// The "default" method is used, beattracking, win_size is the length of FFT
  	// hop_size is the number of frames between two consecutive runs
  	std::string tempo_method = "default";
	if (!(tempo = new_aubio_tempo(const_cast<char *>(tempo_method.c_str()), win_size, hop_size, samplerate))) {
		std::cerr << "new_aubio_tempo failed" << std::endl;
		del_aubio_source(aubio_source);
		aubio_cleanup();
		return false;		
	}

	// Find the tempo
	do {
		// Put data in the input buffer
		aubio_source_do(aubio_source, in, &read);
		// tempo detection
		aubio_tempo_do(tempo, in, out);
		// manipulate the beat
		if (out->data[0] != 0) {
      		std::cout << "beat at " << aubio_tempo_get_last_ms(tempo) << "ms, " << aubio_tempo_get_last_s(tempo) << "s, frame " << aubio_tempo_get_last(tempo) << ", " << aubio_tempo_get_bpm(tempo) << "bpm with confidence " << aubio_tempo_get_confidence(tempo) << std::endl;
      	}
    	n_frames += read;
  	} while (read == hop_size);

  	std::cout << "read " << n_frames * 1. / samplerate << "s, " << n_frames << " frames at " << samplerate << "Hz ("<< n_frames / hop_size << " blocks) from " << file << std::endl;

	// clean up memory
  	del_aubio_tempo(tempo);
  	del_fvec(in);
  	del_fvec(out);
  	del_aubio_source(aubio_source);
  	return true;
}

// http://iub.edu/~emusic/etext/synthesis/chapter4_pv.shtml
bool	getFundamentalFrequencies(const fvec_t *samples) {
	uint_t		window_size = size;
  	fvec_t		*win;
  	fvec_t		*winput;
  	aubio_fft_t	*fft;
  	cvec_t 		*fftOut;

  	if (!(fftOut = new_cvec(window_size))) {
  		std::cerr << "new_cvec for fft failed" << std::endl;
  		return false;
  	}
  	if (!(winput = new_fvec(window_size))) {
  		std::cerr << "new_cvec for winput failed" << std::endl;
  		return false;  		
  	}
  	// Create a envelope called windowing function
  	// The blocks of samples will be multiplied by it
  	// The shape of this window has an effect on the weighting of the resultant analysis
  	// We use the Hanning's method
  	std::string window_method = "hanning";
  	if (!(win = new_aubio_window(const_cast<char *>(window_method.c_str()), window_size))) {
  		std::cerr << "new_aubio_window for fft failed" << std::endl;
  		del_cvec (fftOut);
  		return false;
  	}

  	// Create the fft with the window size
	if (!(fft = new_aubio_fft(window_size))) {
		std::cerr << "new_aubio_fft failed" << std::endl;
		return false;
	}

	// Multiply the block by the envelope (windowing function)
	for (int i = 0; i < window_size; ++i) {
		winput->data[i] = win->data[i] * samples->data[i];
	}


	return true;
}

/**
 * Read the sound file and display a block of hop_size samples
 * @param  file the soung file
 * @return      true if there was no error, otherwise false
 */
bool	readWav(const std::string &file) {
	uint_t 	samplerate = 0;
  	uint_t 	hop_size = 256;
  	uint_t 	n_frames = 0;
  	uint_t	read = 0;

  	aubio_source_t *aubio_source = new_aubio_source(const_cast<char *>(file.c_str()), samplerate, hop_size);
  	if (!aubio_source) {
  		std::cerr << "new_aubio_source failed" << std::endl;
  		return false;
  	}
  	
  	fvec_t *vec = new_fvec(hop_size);
  	samplerate =  aubio_source_get_samplerate(aubio_source);

  	do {
  		aubio_source_do(aubio_source, vec, &read);
  		fvec_print(vec);
  		n_frames += read;
  	} while (read == hop_size);

  	std::cout << "read " << n_frames << " frames at " << samplerate << "Hz (" << n_frames / hop_size << " blocks) from " << file << std::endl;
	del_fvec(vec);
 	del_aubio_source(aubio_source);
  	return true;
}

void	getFeatures(const std::string &file) {
	//readWav(file);
	getBeat(file);
}


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
void	ProcessSoundMidiToWave(std::string const &file, std::string const &directory) {
	std::string extension = GetFileExtension(file);
	if (extension == "mp3") {
		std::string filename = file.substr(0, file.length() - (extension.length() + 1));
		std::string cmd = "timidity -OwM \"" + directory + "/" + file + "\"";
		ExecCommand(cmd);
	}
}


bool	NormalizeDataSet(std::string const &directory) {
	DIR		*dir;
	dirent	*pdir;

	if ((dir = opendir(directory.c_str())) == NULL){
		std::cerr << "Open directory " << directory << " failed" << std::endl;
		return false;
	}

	while ((pdir = readdir(dir))) {
		if (pdir->d_type == DT_REG) {
			 ProcessSoundMidiToWave(pdir->d_name, directory);
		}
	}
	closedir(dir);
	return true;
}

int main(int ac, char **av){
	getFeatures("eiffel65.wav");/*
	if (ac < 2) {
		std::cout << "Usage:" << av[0] << " <training_directory>" << std::endl;
		return 1;
	}
	normalize_dataset(av[1]);*/
    return 0;
}