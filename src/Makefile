CC := gcc
CFLAGS :=  -Wall -Werror -Wextra -std=c11 

GCOVFLAGS := -L. --coverage  -ls21_math
OS = $(shell uname)
ifeq ($(OS), Linux)
	CHECK_FLAGS = $(shell pkg-config --libs check)
else
	CHECK_FLAGS = -lcheck
endif
SRC := $(wildcard s21_*.c)
HDR := $(wildcard s21_*.h)
OBJ := $(patsubst %.c,%.o, ${SRC})

all: clean s21_math.a test

s21_math.o:
		$(CC) $(CFLAGS) -c $(SRC)

s21_math.a: s21_math.o
		ar rcs libs21_math.a $(OBJ)
		ranlib libs21_math.a
		cp libs21_math.a s21_math.a

test: s21_math.a
		$(CC) $(CFLAGS) -o test tests.c $^ $(CHECK_FLAGS)
		./test

gcov_report: test
	$(CC) $(CFLAGS) $(SRC) tests.c libs21_math.a $(CHECK_FLAGS) $(GCOVFLAGS) -o $@
	./gcov_report
	lcov -t gcov_report -o rep.info -c -d .
	genhtml -o report rep.info
	open ./report/index.html
	rm -rf *.gcda *.gcn0 *.info

check:
	cp ./.clang-format .
	clang-format -n *.c *.h
	make test
	CK_FORK=no leaks --atExit -- ./test
	rm .clang-format

clean:
	rm -rf test *.a *.o *.gcno *.gcda *.info test_output gcov_report report
