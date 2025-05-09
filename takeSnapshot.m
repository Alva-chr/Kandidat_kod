function [] = takeSnapshot(suffix, snapshotName)

    filename = snapshotName + "_" + suffix + ".fig";
    fullpath = fullfile(snapshotName, filename);
    saveas(gcf, fullpath);
end