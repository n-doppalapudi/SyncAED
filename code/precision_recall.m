function [precision,recall]=precision_recall(baddata_index,inserted_baddata)
precision=length(intersect(baddata_index,inserted_baddata))/length(baddata_index);
recall=length(intersect(baddata_index,inserted_baddata))/length(inserted_baddata);

end
