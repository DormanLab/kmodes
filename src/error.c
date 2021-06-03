/**
 * @file error.c
 * @author Karin Dorman, kdorman@iastate.edu
 *
 * Functions for outputting error and debugging messages.
 */

#include "error.h"

int global_debug_level = ABSOLUTE_SILENCE;	/* allow some informational */

int vmessage(FILE *fp, const char *file_name, const char *fxn_name, int line,
	int msg_type, int msg_id, const char *msg, va_list vl);

/**
 * Write uniformly formatted message to requested stream.
 * Write uniformly formatted debug or error messages to requested stream.
 * Certain common errors/problems issue default messages so they need not be
 * typed repeatedly and nonstandardly.  Returns msg_id so that the function
 * can be conveniently called as return message(...) or exit(message(...)).
 *
 * @param fp file handle where message to be written
 * @param file_name name of offending source file
 * @param fxn_name name of offending function
 * @param line line number of offending function
 * @param msg_type urgency of message
 * @param msg_id identify one of a default set of canned messages
 * @param msg formatted string to output via vprintf
 * @return msg_id
 */
int message(FILE *fp, const char *file_name, const char *fxn_name, int line,
	int msg_type, int msg_id, char const * const msg, ...) {
	int ret;
	va_list vl;
	va_start(vl, msg);
	ret = vmessage(fp, file_name, fxn_name, line, msg_type, msg_id, msg, vl);
	va_end(vl);
	return(ret);
} /* message */

int vmessage(FILE *fp, const char *file_name, const char *fxn_name, int line,
	int msg_type, int msg_id, const char *msg, va_list vl) {
	int nsec;
	kmodes_fprintf(fp, "%s [%s::%s(%4d)]: ",
		msg_type == INFO_MSG ? "INFO" : msg_type == DEBUG_MSG ? "DEBUG"
			: msg_type == WARNING_MSG ? "WARNING" : "ERROR",
		file_name, fxn_name, line);
	if (msg_id == NO_ERROR) {
		kmodes_vfprintf(fp, msg, vl);
	} else {
		switch(msg_id) {
			case MEMORY_ALLOCATION:
				if (msg) {
					kmodes_fprintf(fp, "could not allocate ");
					kmodes_vfprintf(fp, msg, vl);
				} else
					kmodes_fprintf(fp, "memory allocation error\n");
				break;
			case INVALID_CMD_OPTION:
				kmodes_fprintf(fp, "unrecognized command option");
				if (msg) {
					kmodes_fprintf(fp, ": ");
					kmodes_vfprintf(fp, msg, vl);
				} else
					kmodes_fprintf(fp, "\n");
				break;
			case INVALID_CMD_ARGUMENT:
				kmodes_fprintf(fp, "invalid argument to command option");
				if (msg) {
					kmodes_fprintf(fp, ": ");
					kmodes_vfprintf(fp, msg, vl);
				} else
					kmodes_fprintf(fp, "\n");
				break;
			case INVALID_CMDLINE:
				kmodes_fprintf(fp, "[invalid command line] ");
				if (msg)
					kmodes_vfprintf(fp, msg, vl);
				else
					kmodes_fprintf(fp, "\n");
				break;
			case INVALID_USER_INPUT:
				kmodes_fprintf(fp, "[invalid user choice] ");
				if (msg)
					kmodes_vfprintf(fp, msg, vl);
				else
					kmodes_fprintf(fp, "\n");
				break;
			case FILE_OPEN_ERROR:
				kmodes_fprintf(fp, "could not open file \"%s\"\n", msg);
				break;
			case PIPE_OPEN_ERROR:
				kmodes_fprintf(fp, "could not execute command \"%s\"\n", msg);
				break;
			case FILE_NOT_FOUND:
				kmodes_fprintf(fp, "file \"%s\" not found\n", msg);
				break;
			case FILE_FORMAT_ERROR:
				kmodes_fprintf(fp, "invalid file format");
				if (msg) {
					kmodes_fprintf(fp, ": ");
					kmodes_vfprintf(fp, msg, vl);
				} else
					kmodes_fprintf(fp, "\n");
				break;
			case END_OF_FILE:
				kmodes_fprintf(fp, "unexpected end of file in file \"%s\"\n", msg);
				break;
			case INTERNAL_MISMATCH:
				kmodes_fprintf(fp, "[internal mismatch] %s\n", msg);
				break;
			case OUT_OF_TIME:
				kmodes_fprintf(fp, "out of time");
				if (msg) {
					nsec = va_arg(vl, int);
					kmodes_fprintf(fp, " (%slimit %02d:%02dm)",
						msg, (int)(nsec/3600),
						(int)((nsec%3600)/60));
				}
				kmodes_fprintf(fp, "\n");
				break;
			case MEMORY_USAGE_LIMIT:
				kmodes_fprintf(fp, "exceed memory limit");
				if (msg) {
					kmodes_fprintf(fp, ": ");
					kmodes_vfprintf(fp, msg, vl);
				} else
					kmodes_fprintf(fp, "\n");
				break;
			case EXCEED_ITERATIONS:
				kmodes_fprintf(fp, "exceed iteration limit");
				if (msg) {
					kmodes_fprintf(fp, ": ");
					kmodes_vfprintf(fp, msg, vl);
				} else
					kmodes_fprintf(fp, "\n");
				break;
			default:
				kmodes_vfprintf(fp, msg, vl);
		}
	}
	return msg_id;
} /* message */

#ifdef USE_CURSES
int wmessage(WINDOW *wp, const char *file_name, const char *fxn_name, int line,
	int msg_type, int msg_id, const char *msg, ...) {
	int nsec;
	va_list args;
	wprintw(wp, "%s [%s::%s(%d)]: ",
		msg_type == INFO_MSG ? "INFO" : msg_type == DEBUG_MSG ? "DEBUG"
			: msg_type == WARNING_MSG ? "WARNING" : "ERROR",
		file_name, fxn_name, line);
	if (msg_id == NO_ERROR) {
		va_start(args, msg);
		vw_printw(wp, msg, args);
		va_end(args);
	} else {
		switch(msg_id) {
			case MEMORY_ALLOCATION:
				if (msg) {
					wprintw(wp, "could not allocate ");
					va_start(args, msg);
					vw_printw(wp, msg, args);
					va_end(args);
				} else
					wprintw(wp, "memory allocation error\n");
				break;
			case INVALID_CMD_OPTION:
				wprintw(wp, "unrecognized command option");
				if (msg) {
					wprintw(wp, ": ");
					va_start(args, msg);
					vw_printw(wp, msg, args);
					va_end(args);
				} else
					wprintw(wp, "\n");
				break;
			case INVALID_CMD_ARGUMENT:
				wprintw(wp, "invalid argument to command option");
				if (msg) {
					wprintw(wp, ": ");
					va_start(args, msg);
					vw_printw(wp, msg, args);
					va_end(args);
				} else
					wprintw(wp, "\n");
				break;
			case INVALID_CMDLINE:
				wprintw(wp, "[invalid command line] %s\n", msg);
				break;
			case INVALID_USER_INPUT:
				wprintw(wp, "[invalid user choice] %s\n", msg);
				break;
			case FILE_OPEN_ERROR:
				wprintw(wp, "could not open file \"%s\"\n", msg);
				break;
			case FILE_NOT_FOUND:
				wprintw(wp, "file \"%s\" not found\n", msg);
				break;
			case FILE_FORMAT_ERROR:
				wprintw(wp, "invalid file format");
				if (msg) {
					wprintw(wp, ": ");
					va_start(args, msg);
					vw_printw(wp, msg, args);
					va_end(args);
				} else
					wprintw(wp, "\n");
				break;
			case END_OF_FILE:
				wprintw(wp, "unexpected end of file in file \"%s\"\n", msg);
				break;
			case INTERNAL_MISMATCH:
				wprintw(wp, "[internal mismatch] %s\n", msg);
				break;
			case OUT_OF_TIME:
				wprintw(wp, "out of time");
				if (msg) {
					va_start(args, msg);
					nsec = va_arg(args, int);
					va_end(args);
					wprintw(wp, " (%slimit %02d:%02dm)",
						msg, (int)(nsec/3600),
						(int)((nsec%3600)/60));
				}
				wprintw(wp, "\n");
				break;
			case MEMORY_USAGE_LIMIT:
				wprintw(wp, "exceed memory limit");
				if (msg) {
					wprintw(wp, ": ");
					va_start(args, msg);
					vw_printw(wp, msg, args);
					va_end(args);
				} else
					wprintw(wp, "\n");
				break;
			case EXCEED_ITERATIONS:
				wprintw(wp, "exceed iteration limit");
				if (msg) {
					wprintw(wp, ": ");
					va_start(args, msg);
					vw_printw(wp, msg, args);
					va_end(args);
				} else
					wprintw(wp, "\n");
				break;
			default:
				va_start(args, msg);
				vw_printw(wp, msg, args);
				va_end(args);
		}
	}
	return msg_id;
} /* wmessage */
#endif
