licenses(["notice"])

cc_library(
    name = "Dendrogram",
    hdrs = ["Dendrogram.h"],
    deps = [
        ":rctree_hashtable",
        ":rctree_utils",
        "@gbbs//benchmarks/Connectivity/SimpleUnionAsync:Connectivity",
        "@gbbs//gbbs",
        "@gbbs//gbbs:julienne",
    ],
)

cc_library(
    name = "concurrent_table",
    hdrs = ["concurrent_table.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "rctree_utils",
    hdrs = ["rctree_utils.h"],
    deps = [
        "@gbbs//gbbs",
    ],
)

cc_library(
    name = "rctree_hashtable",
    hdrs = ["rctree_hashtable.h"],
    deps = [
        ":concurrent_table",
        ":rctree_utils",
        "@gbbs//gbbs",
    ],
)

cc_binary(
    name = "Dendrogram_main",
    srcs = ["Dendrogram.cc"],
    deps = [
        ":Dendrogram",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)