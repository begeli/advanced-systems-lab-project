#pragma once

/// Provides useful things for testing
#include <cstdio>
#include <iostream>
#include <gtest/gtest.h>
#include <fcntl.h>

/// Capture std::cout output
class CaptureCout {
 public:
  std::streambuf* sbuf;
  std::stringstream stream;

  CaptureCout() : sbuf(std::cout.rdbuf()) {
    std::cout.rdbuf(stream.rdbuf());
  }

  ~CaptureCout() {
    std::cout.rdbuf(sbuf);
  }
};

// source: https://stackoverflow.com/questions/4832603/how-could-i-temporary-redirect-stdout-to-a-file-in-a-c-program/4832902#4832902
class CaptureStdout {
 private:
  int stdout_fd;
 public:
  // TODO figure out how to read from the file we redirect to
  // then we can compare printf output with an expected string
  // Replace /dev/null with whatever file you want to redirect stdout to.
  CaptureStdout() :stdout_fd(dup(STDOUT_FILENO)) {
    int redir_fd = open("/dev/null", O_WRONLY);
    dup2(redir_fd, STDOUT_FILENO);
    close(redir_fd);
  }

  ~CaptureStdout() {
    fflush(stdout);
    dup2(stdout_fd, STDOUT_FILENO);
    close(stdout_fd);
  }
};
