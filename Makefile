CC?=cc
CFLAGS= -O3 -Wall -Wpedantic -lm -mavx
NAME=spinny
SRC=spinny.c
SDL_FLAGS!=sdl2-config --cflags --libs 

${NAME}: ${SRC}
	${CC} ${CFLAGS} ${SDL_FLAGS}  ${SRC} -o ${NAME}

clean:
	rm ${NAME}

.PHONY: clean
