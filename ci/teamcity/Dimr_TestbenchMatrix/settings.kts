import jetbrains.buildServer.configs.kotlin.*

import testbenchMatrix.Trigger
import testbenchMatrix.Linux
import testbenchMatrix.Windows
import testbenchMatrix.Release

version = "2024.03"

project {
    description = "contact: BlackOps (black-ops@deltares.nl)"

    buildType(Trigger)
    buildType(Linux)
    buildType(Windows)
    buildType(Release)

    buildTypesOrder = arrayListOf(Trigger, Linux, Windows, Release)
}