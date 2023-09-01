licenses(["notice"])

cc_library(
    name = "Dendrogram",
    hdrs = ["Dendrogram.h"],
    deps = [
        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "@gbbs//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "@gbbs//gbbs",
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
