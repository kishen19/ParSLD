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
    name = "caterpillar",
    hdrs = ["caterpillar.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "star",
    hdrs = ["star.h"],
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
    name = "graph_io",
    hdrs = ["graph_io.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "generate_MSF",
    hdrs = ["generate_MSF.h"],
    deps = [
        ":graph_io",
        "@gbbs//gbbs",
        "@gbbs//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "@gbbs//benchmarks/MinimumSpanningForest/Kruskal:MinimumSpanningForest",
    ],
)

cc_library(
    name = "benchmark",
    hdrs = ["benchmark.h"],
    deps = [
        "@gbbs//gbbs",
        "@gbbs//benchmarks/TriangleCounting/ShunTangwongsan15:Triangle",
        "@gbbs//gbbs:graph_io",
        ":caterpillar",
        ":generate_MSF",
        ":graph_io",
        ":path",
        ":star",
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
