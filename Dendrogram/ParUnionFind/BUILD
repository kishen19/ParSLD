licenses(["notice"])

cc_library(
    name = "Dendrogram",
    hdrs = ["Dendrogram.h"],
    deps = [
        "//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "//gbbs",
    ],
)

cc_binary(
    name = "Dendrogram_main",
    srcs = ["Dendrogram.cc"],
    deps = [":Dendrogram"],
)

package(
    default_visibility = ["//visibility:public"],
)
