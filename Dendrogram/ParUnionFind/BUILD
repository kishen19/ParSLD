licenses(["notice"])

cc_library(
    name = "DendAsync",
    hdrs = ["DendAsync.h"],
    deps = [
        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "@gbbs//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "@gbbs//gbbs",
        "//Dendrogram/heaps:leftist_heap",
        "//Dendrogram/heaps:skew_heap",
        "//Dendrogram/heaps:pairing_heap",
    ],
)

cc_library(
    name = "DendSync",
    hdrs = ["DendSync.h"],
    deps = [
        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "@gbbs//benchmarks/MinimumSpanningForest/Boruvka:MinimumSpanningForest",
        "@gbbs//gbbs",
        "//Dendrogram/heaps:leftist_heap",
        "//Dendrogram/heaps:skew_heap",
        "//Dendrogram/heaps:pairing_heap",
    ],
)

cc_binary(
    name = "Dendrogram_main",
    srcs = ["Dendrogram.cc"],
    deps = [
        ":DendAsync",
        ":DendSync",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
