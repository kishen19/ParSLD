licenses(["notice"])

cc_library(
    name = "full_binary",
    hdrs = ["full_binary.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "path",
    hdrs = ["path.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "uniform_hook",
    hdrs = ["uniform_hook.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "benchmark",
    hdrs = ["benchmark.h"],
    deps = [
        "@gbbs//gbbs",
        ":path",
        ":full_binary",
        ":uniform_hook",
    ],
)

cc_binary(
    name = "BuildFullBinaryTree",
    srcs = ["full_binary.cc"],
    deps = [
        "@gbbs//gbbs",
        ":full_binary",
    ],
)

cc_binary(
    name = "BuildPath",
    srcs = ["path.cc"],
    deps = [
        "@gbbs//gbbs",
        ":path",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
