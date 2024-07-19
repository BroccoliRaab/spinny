CC?=cc
CFLAGS=-g -O0 -Wall -Wpedantic
NAME=spinny
SRC=spinny.c
SDL_FLAGS!=sdl2-config --cflags --libs --static-libs

${NAME}: ${SRC}
	${CC} ${CFLAGS} ${SDL_FLAGS}  ${SRC} -o ${NAME}

clean:
	rm ${NAME}

.PHONY: clean
