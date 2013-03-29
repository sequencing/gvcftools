#
# IDS must be defined:
#
IDS := CXXFLAGS LDFLAGS LDLIBS LIBUTIL

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

TEST_SRCS = $(wildcard test/*.cpp)
TEST_OBJS = $(TEST_SRCS:%.cpp=%.o)

.PHONY: clean test


$(LIBUTIL): $(OBJS)
	$(AR) -csru $@ $(OBJS)

test: test_execute

TEST_PROGRAM := $(LIBUTIL)_unittest

test_execute: $(TEST_PROGRAM)
	@echo Unit testing $(LIBUTIL)
	@$(TEST_PROGRAM) --log_level=test_suite
	#@$(TEST_PROGRAM) --log_level=all

#$(TEST_PROGRAM): LDLIBS += -lboost_unit_test_framework

$(TEST_PROGRAM): $(LIBUTIL) $(TEST_OBJS)
	$(CXX) $(TEST_OBJS) $(LIBUTIL) -o $(TEST_PROGRAM) $(LDFLAGS) $(LDLIBS) -lboost_unit_test_framework

$(TEST_OBJS): CXXFLAGS += -I$(CURDIR)

clean:
	rm -f $(LIBUTIL) $(OBJS) $(TESTPROG) $(TEST_OBJS)

