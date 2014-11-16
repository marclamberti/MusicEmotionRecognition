CC	= g++ -std=c++11
CFLAGS 	= -c -Wall
SRCS	= $(wildcard *.cc)
OBJS	= $(SRCS:.cc=.o)
BIN	= merge_results

all: $(SRCS) $(BIN)

$(BIN): $(OBJS)
		$(CC) $(LDFLAGS) $(OBJS) -o $@

.cc.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
		rm -rf $(OBJS)

fclean: clean
		rm -rf $(BIN)

re: fclean all