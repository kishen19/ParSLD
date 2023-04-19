licenses(["notice"])

cc_library(
    name = "Dendogram",
    hdrs = ["Dendogram.h"],
    deps = [
        "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "//gbbs",
    ],
)

cc_binary(
    name = "Dendogram_main",
    srcs = ["Dendogram.cc"],
    deps = [":Dendogram"],
)

package(
    default_visibility = ["//visibility:public"],
)
