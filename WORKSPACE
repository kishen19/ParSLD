load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "new_git_repository")
load("@bazel_tools//tools/cpp:cc_configure.bzl", "cc_configure")

cc_configure()

local_repository(
   name = "gbbs",
   path = "external/gbbs",
)

local_repository(
   name = "fast-bcc",
   path = "external/FAST-BCC",
)

http_archive(
    name = "parlaylib",
    sha256 = "68c062ad116fd49d77651d7a24fb985aa66e8ec9ad05176b6af3ab5d29a16b1f",
    strip_prefix = "parlaylib-bazel/include/",
    urls = ["https://github.com/ParAlg/parlaylib/archive/refs/tags/bazel.tar.gz"],
)

http_archive(
    name = "googletest",
    sha256 = "b4870bf121ff7795ba20d20bcdd8627b8e088f2d1dab299a031c1034eddc93d5",
    strip_prefix = "googletest-release-1.11.0",
    urls = ["https://github.com/google/googletest/archive/release-1.11.0.tar.gz"],
)